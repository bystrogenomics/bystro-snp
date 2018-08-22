package main

import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"runtime/pprof"
	"strconv"
	"strings"
	"sync"

	"github.com/akotlar/bystro-utils/parse"
)

const (
	concurrency int = 3
	chromIdx int = 0
	posIdx int = 1
	refIdx int = 2
	altIdx int = 3
	altCountIdx int = 4
	typeIdx int = 5
	firstSampleIdx int = 6
)

const (
	tabByte = byte('\t')
	clByte = byte('\n')
	zeroByte = byte('0')
)

var fileMutex sync.Mutex

// Arrays of length 2 are hets, of length 1 are homozygote
var iupac = map[byte][]byte{
	'A': []byte{'A'}, 'C': []byte{'C'}, 'G': []byte{'G'}, 'T': []byte{'T'}, 'D': []byte{'-'}, 'I': []byte{'+'},
	'R': []byte{'A', 'G'}, 'Y': []byte{'C', 'T'}, 'S': []byte{'G', 'C'},
	'W': []byte{'A', 'T'}, 'K': []byte{'G', 'T'}, 'M': []byte{'A', 'C'},
	// Fake; there is no other allele
	'E': []byte{'-'}, 'H': []byte{'+'},
}

type Config struct {
	inPath         string
	outPath        string
	errPath        string
	cpuProfile     string
	emptyField     string
	fieldDelimiter string
	minGq          float64
}

func setup(args []string) *Config {
	config := &Config{}
	flag.StringVar(&config.inPath, "in", "", "The input file path (optional: default is stdin)")
	flag.StringVar(&config.outPath, "out", "", "The output file path (optional: default is stdout)")
	flag.StringVar(&config.cpuProfile, "cProf", "", "Path to output cpu profile")
	flag.StringVar(&config.errPath, "errPath", "", "The output path for the JSON output (optional)")
	flag.StringVar(&config.emptyField, "emptyField", "!", "The output path for the JSON output (optional)")
	flag.StringVar(&config.fieldDelimiter, "fieldDelimiter", ";", "The output path for the JSON output (optional)")
	flag.Float64Var(&config.minGq, "minGq", .95, "The minimum confidence of a genotype call")

	// allows args to be mocked https://github.com/nwjlyons/email/blob/master/inputs.go
	// can only run 1 such test, else, redefined flags error
	a := os.Args[1:]
	if args != nil {
		a = args
	}

	flag.CommandLine.Parse(a)

	return config
}

func init() {
	log.SetFlags(0)
}

func main() {
	config := setup(nil)

	inFh := (*os.File)(nil)
	if config.inPath != "" {
		var err error

		inFh, err = os.Open(config.inPath)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		inFh = os.Stdin
	}

	// make sure it gets closed
	defer inFh.Close()
	if config.errPath != "" {
		var err error
		os.Stderr, err = os.Open(config.errPath)
		if err != nil {
			log.Fatal(err)
		}
	}

	outFh := (*os.File)(nil)
	if config.outPath != "" {
		var err error

		outFh, err = os.OpenFile(config.outPath, os.O_WRONLY|os.O_CREATE, 0664)
		if err != nil {
			log.Fatal(err)
		}
	} else {
		outFh = os.Stdout
	}
	// make sure it gets closed
	defer outFh.Close()

	if config.cpuProfile != "" {
		f, err := os.Create(config.cpuProfile)
		if err != nil {
			log.Fatal(err)
		}
		pprof.StartCPUProfile(f)
		defer pprof.StopCPUProfile()
	}

	reader := bufio.NewReader(inFh)

	writer := bufio.NewWriter(outFh)

	fmt.Fprintln(writer, strings.Join(parse.Header, "\t"))

	readSnp(config, reader, writer)

	writer.Flush()
	outFh.Close()
}

func readSnp(config *Config, reader *bufio.Reader, writer *bufio.Writer) {
	// Read buffer
	workQueue := make(chan string, 100)
	complete := make(chan bool)

	endOfLineByte, numChars, headerLine, err := parse.FindEndOfLine(reader, "")

	if err != nil {
		log.Fatal(err)
	}

	header := strings.Split(headerLine[:len(headerLine)-numChars], "\t")

	if len(header) == 0 {
		log.Fatal("No header found")
	}

	// Remove periods from sample names
	parse.NormalizeHeader(header)

	// Read the lines into the work queue.
	go func() {
		for {
			row, err := reader.ReadString(endOfLineByte) // 0x0A separator = newline

			if err == io.EOF {
				break
			} else if err != nil {
				log.Fatal(err)
			} else if row == "" {
				// We may have not closed the pipe, but not have any more information to send
				// Wait for EOF
				continue
			}

			workQueue <- row[:len(row)-numChars]
		}

		// Close the channel so everyone reading from it knows we're done.
		close(workQueue)
	}()

	// Now read them all off, concurrently.
	for i := 0; i < concurrency; i++ {
		go processLine(header, config, workQueue, writer, complete)
	}

	// Wait for everyone to finish.
	for i := 0; i < concurrency; i++ {
		<-complete
	}

}

func processLine(header []string, config *Config,
	queue chan string, writer *bufio.Writer, complete chan bool) {

	emptyField := config.emptyField
	fieldDelimiter := config.fieldDelimiter
	minGq := config.minGq

	alleleCache := make(map[byte]map[string][]string)
	var altAlleles []string

	var homs [][]string
	var hets [][]string
	var missing [][]string
	var effectiveSamples float64

	var numSamples float64

	//"Header is not even, last genotype placeholder was chopped, adding back 1 field"
	if (len(header)-firstSampleIdx)%2 != 0 {
		numSamples = float64(len(header)+1-firstSampleIdx) / 2
	} else {
		numSamples = float64(len(header)-firstSampleIdx) / 2
	}

	if numSamples != float64(int64(numSamples)) {
		log.Printf("Warning: calculated # samples is %f, but expected to be integer", numSamples)
	}

	var output bytes.Buffer
	var oIdx int
	for line := range queue {
		record := strings.Split(line, "\t")

		if !validType(record[typeIdx]) {
			continue
		}

		if oIdx >= 10000 {
			fileMutex.Lock()

			writer.Write(output.Bytes())

			fileMutex.Unlock()

			output.Reset()

			oIdx = 0
		}

		// Hit rate for this is quite high; runs only 12 times for all SNPs
		// and N times for indels, multiallelics, but those are rare, and do repeat
		altAlleles = gatherAlt(record[refIdx][0], record[altIdx], alleleCache)

		if numSamples > 0 {
			homs, hets, missing = makeHetHomozygotes(record, header, altAlleles, minGq)
		}

		for i, alt := range altAlleles {
			// If no samples are provided, annotate what we can, skipping hets and homs
			// If we have samples, but no non-missing found at this site, skip the site
			if numSamples > 0 {
				// this site has no samples at all with the minor allele, so skip it
				if len(hets[i]) == 0 && len(homs[i]) == 0 {
					continue
				}

				// Missing values should all be the same; currently makeHetHomozygotes
				// provides array of N (equal) missingngess arrays for convenience, one for each alt
				// to match the use of hets and homs
				// This use will incur minor hit for multiallelics, but leaves open
				// possibility that missing[0] != missing[1] at some point in future
				if len(missing) > 0 && len(missing[i]) > 0 {
					effectiveSamples = numSamples - float64(len(missing[i]))
				} else {
					effectiveSamples = numSamples
				}
			}

			output.WriteString(record[chromIdx])
			output.WriteByte(tabByte)
			output.WriteString(record[posIdx])
			output.WriteByte(tabByte)

			if len(altAlleles) > 1 {
				output.WriteString("MULTIALLELIC")
			} else {
				output.WriteString(record[typeIdx])
			}

			output.WriteByte(tabByte)
			output.WriteString(record[refIdx])
			output.WriteByte(tabByte)
			output.WriteString(alt)
			output.WriteByte(tabByte)

			if len(altAlleles) > 1 {
				output.WriteString(parse.NotTrTv)
			} else {
				output.WriteString(parse.GetTrTv(record[refIdx], alt))
			}

			output.WriteByte(tabByte)

			if len(hets) == 0 || len(hets[i]) == 0 {
				output.WriteString(emptyField)
				output.WriteByte(tabByte)
				output.WriteByte(zeroByte)
			} else {
				output.WriteString(strings.Join(hets[i], fieldDelimiter))
				output.WriteByte(tabByte)

				// This gives plenty precision; we are mostly interested in
				// the first or maybe 2-3 significant digits
				// https://play.golang.org/p/Ux-QmClaJG
				// Also, gnomAD seems to use 6 bits of precision
				// the bitSize == 64 allows us to round properly past 6 s.f
				// Note: 'G' requires these numbers to be < 0 for proper precision
				// (elase only 6 s.f total, rather than after decimal)

				// heterozygosity is relative to the number of complete samples
				output.WriteString(strconv.FormatFloat(float64(len(hets[i]))/effectiveSamples, 'G', 4, 64))
			}

			output.WriteByte(tabByte)

			if len(homs) == 0 || len(homs[i]) == 0 {
				output.WriteString(emptyField)
				output.WriteByte(tabByte)
				output.WriteByte(zeroByte)
			} else {
				output.WriteString(strings.Join(homs[i], fieldDelimiter))
				output.WriteByte(tabByte)

				// homozygosity is relative to the number of complete samples
				output.WriteString(strconv.FormatFloat(float64(len(homs[i]))/effectiveSamples, 'G', 4, 64))
			}

			output.WriteByte(tabByte)

			// Missing values should all be the same; currently makeHetHomozygotes
			// provides array of N (equal) missingngess arrays for convenience, one for each alt
			// to match the use of hets and homs
			if len(missing) == 0 || len(missing[i]) == 0 {
				output.WriteString(emptyField)
				output.WriteByte(tabByte)
				output.WriteByte(zeroByte)
			} else {
				output.WriteString(strings.Join(missing[i], fieldDelimiter))
				output.WriteByte(tabByte)

				// missingness is relative to the total number of samples
				output.WriteString(strconv.FormatFloat(float64(len(missing[i]))/numSamples, 'G', 4, 64))
			}

			// Write the sample minor allele frequency
			// This can be 0 in one of wo situations
			// First, if we have only missing genotypes at this site
			// However, in this case, we don't reach this code, because of line
			// 302 (if len(homs) == 0 && len(hets) == 0)
			// Else if there are truly no minor allele
			output.WriteByte(tabByte)

			//For sites without samples
			if effectiveSamples == 0 {
				output.WriteByte(zeroByte)
			} else {
				// sampleMaf numerator is just 2x the number of homozygotes + number of heterozygotes
				// because .snp files list haploids exactly the same as homozygous diploids
				// sampleMaf denominator is typically (len(fields) - len(missing))*2 - 6
				// where 6 is the number of fields before the first sample field (or equivalently
				// it's the first sample field index)
				// Note that len(fields) - 6 contains 2x as many fields as samples, with every other
				// being a genotype
				// This nicely matches our expectation for the number of alleles per sample
				// allowing us to just subtract len(missing) * 2
				// Also note that we use len(fields) not len(header)
				// Depending on the reader implementation, header may contain 1 fewer than the expected
				// number of fields, because the last character is an empty string, followed by a "\n"
				// So, for instance, Perl's chomp will remove not only the "\n", but also rewind the record
				// the field on the left side of the preceeding "\t"
				output.WriteString(strconv.FormatFloat(float64(len(homs[i])*2+len(hets[i]))/(effectiveSamples*2), 'G', 4, 64))
			}

			output.WriteByte(clByte)
		}
	}

	if output.Len() > 0 {
		fileMutex.Lock()

		writer.Write(output.Bytes())

		fileMutex.Unlock()
	}

	// log.Println("Worker hit, missed this many times: ", hitCount, missCount)
	// Let the main process know we're done.
	complete <- true
}

// Explicitly whitelist call types, because this gets stored, could be malicious
func validType(cType string) bool {
	return cType == "SNP" || cType == "INS" || cType == "DEL" || cType == "MULTIALLELIC" || cType == "DENOVO_SNP" || cType == "DENOVO_INS" || cType == "DENOVO_DEL" || cType == "DENOVO_MULTIALLELIC"
}

// Save effort in identifying alternate alleles
func gatherAlt(ref byte, alleles string, alt map[byte]map[string][]string) []string {
	if alt[ref] == nil {
		alt[ref] = make(map[string][]string)
	} else if alt[ref][alleles] != nil {
		return alt[ref][alleles]
	}

	if !strings.Contains(alleles, ",") {
		alt[ref][alleles] = []string{alleles}
		return alt[ref][alleles]
	}

	for _, val := range strings.Split(alleles, ",") {
		if val[0] == ref {
			continue
		}

		alt[ref][alleles] = append(alt[ref][alleles], val)
	}

	return alt[ref][alleles]
}

// TODO: Right now we define missing to always be the same across all alleles in a multiallelic
// since such sites are necessarily ambiguous, we can't assign them to a particular allele
// However, this funciton defines missing to be 2D array; we can optimize this away
func makeHetHomozygotes(fields []string, header []string, altAlleles []string, minGQ float64) ([][]string, [][]string, [][]string) {
	missing := make([][]string, len(altAlleles), len(altAlleles))
	het := make([][]string, len(altAlleles), len(altAlleles))
	hom := make([][]string, len(altAlleles), len(altAlleles))

SAMPLES:
	for i := firstSampleIdx; i < len(fields); i += 2 {
		if fields[refIdx] == fields[i] {
			continue
		}

		if fields[i] == "N" {
			parse.AppendMissing(len(altAlleles), header[i], missing)
			continue
		}

		conf, err := strconv.ParseFloat(fields[i+1], 64)

		if err != nil {
			log.Printf("%s:%s: %s genotype invalid confidence %s", fields[chromIdx], fields[posIdx], header[i], fields[i+1])
			parse.AppendMissing(len(altAlleles), header[i], missing)
			continue
		}

		if conf < minGQ {
			parse.AppendMissing(len(altAlleles), header[i], missing)
			continue
		}

		//Fast path for homozygotes
		//avoid calculating hash; may be faster for small sets
		//small hash linear search golang
		//GC somewhat depleted in human http://blog.kokocinski.net/index.php/gc-content-of-human-chromosomes?blog=2
		if fields[i] == "A" || fields[i] == "T" || fields[i] == "C" || fields[i] == "G" || fields[i] == "D" || fields[i] == "I" {
			for altIndex, oAlt := range altAlleles {
				if fields[i] == oAlt || (oAlt[0] == '-' && fields[i] == "D") || (oAlt[0] == '+' && fields[i] == "I") {
					hom[altIndex] = append(hom[altIndex], header[i])
					continue SAMPLES
				}
			}

			parse.AppendMissing(len(altAlleles), header[i], missing)
			log.Printf("%s:%s: %s genotype %s not in Alleles", fields[chromIdx], fields[posIdx], header[i], fields[i])
			continue
		}

		// Heterozygote
		iupacArr, ok := iupac[fields[i][0]]

		if !ok {
			log.Printf("%s:%s: %s genotype %s not IUPAC", fields[chromIdx], fields[posIdx], header[i], fields[i])
			parse.AppendMissing(len(altAlleles), header[i], missing)
			continue
		}

	IUPAC:
		for _, tAlt := range iupacArr {
			if tAlt == fields[refIdx][0] {
				continue
			}

			// First check that the iupac code makes sense given the present alleles
			// It could be that (in low-enough quality sites) the code doesn't actually correspond
			// to the present alleles
			for altIndex, oAlt := range altAlleles {
				if tAlt == oAlt[0] {
					het[altIndex] = append(het[altIndex], header[i])
					continue IUPAC
				}
			}

			// If we're here, a sample has a genotype that isn't in Alleles
			// Since genotype error add it to missingGenos, remove from het list if present
			for altIndex, _ := range altAlleles {
				// If previously added to hets, it would be the last item
				if len(het[altIndex]) > 0 && het[altIndex][len(het[altIndex])-1] == header[i] {
					het[altIndex] = het[altIndex][:len(het[altIndex])-1]
				}

				missing[altIndex] = append(missing[altIndex], header[i])
			}

			log.Printf("%s:%s: %s genotype %s not in Alleles", fields[chromIdx], fields[posIdx], header[i], fields[i])
			break IUPAC
		}
	}

	return hom, het, missing
}

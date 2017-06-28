// estab exports elasticsearch fields as tab separated values
package main
import (
	"bufio"
	"bytes"
	"flag"
	"fmt"
	"io"
	"log"
	"os"
	"strings"
	"sync"
	"regexp"
	"strconv"
)
const chromIdx int = 0
const posIdx int = 1
const refIdx int = 2
const altIdx int = 3
const altCountIdx int = 4
const typeIdx int = 5
const firstSampleIdx int = 6

// Arrays of length 2 are hets, of length 1 are homozygote
var iupac = map[byte][]byte {
	'A': []byte{'A'}, 'C': []byte{'C'}, 'G': []byte{'G'}, 'T': []byte{'T'}, 'D': []byte{'-'}, 'I': []byte{'+'},
	'R': []byte{'A','G'}, 'Y': []byte{'C', 'T'}, 'S': []byte{'G', 'C'},
	'W': []byte{'A', 'T'}, 'K': []byte{'G', 'T'}, 'M': []byte{'A', 'C'},
	// Fake; there is no other allele
	'E': []byte{'-'}, 'H': []byte{'+'},
}

type Config struct {
	inPath string
	errPath string
	emptyField string
	fieldDelimiter string
	minGq float64
}

func setup(args []string) *Config {
	config := &Config{}
	flag.StringVar(&config.inPath, "inPath", "", "The input file path (optional: default is stdin)")
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

// NOTE: For now this only supports \n end of line characters
// If we want to support CLRF or whatever, use either csv package, or set a different delimiter
func main() {
	config := setup(nil)

	readSNP(config)
}

func readSNP (config *Config) {
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

	reader := bufio.NewReader(inFh)

	endOfLineByte, numChars, headerLine, err := findEndOfLineChar(reader, "")

	if err != nil {
		log.Fatal(err)
	}

	header := strings.Split(headerLine[:len(headerLine) - numChars], "\t")

	if len(header) == 0 {
		log.Fatal("No header found")
	}

	// Remove periods from sample names
	normalizeSampleNames(header)

	c := make(chan string)

	// I think we need a wait group, not sure.
	wg := new(sync.WaitGroup)

	go func() {
		var record []string
		var altAlleles []byte

		altCache := make(map[string]map[string][]byte)

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

			// remove the trailing \n or \r
			// equivalent of chomp https://groups.google.com/forum/#!topic/golang-nuts/smFU8TytFr4
			record = strings.Split(row[:len(row) - numChars], "\t")

			if !validType(record[typeIdx]) {
				continue
			}

			altAlleles = gatherAlt(record[refIdx], record[altIdx], altCache)

			if len(altAlleles) == 0 {
				continue
			}

			wg.Add(1)
			go processLine(record, header, altAlleles, config.emptyField,
				config.fieldDelimiter, config.minGq, c, wg)
		}

		wg.Wait()
		close(c)
	}()

	// Write all the data
	// Somewhat surprisingly this is faster than building up array and writing in builk
	for data := range c {
		fmt.Print(data)
	}
}

func findEndOfLineChar (r *bufio.Reader, s string) (byte, int, string, error) {
	runeChar, _, err := r.ReadRune()

	if err != nil {
		return byte(0), 0, "", err
	}

	if runeChar == '\r' {
		nextByte, err := r.Peek(1)

		if err != nil {
			return byte(0), 0, "", err
		}

		if rune(nextByte[0]) == '\n' {
			//Remove the line feed
			_, _, err = r.ReadRune()

			if err != nil {
				return byte(0), 0, "", err
			}

			return nextByte[0], 2, s, nil
		}

		return byte('\r'), 1, s, nil
	}

	if runeChar == '\n' {
		return byte('\n'), 1, s, nil
	}

	s += string(runeChar)
	return findEndOfLineChar(r, s)
}

// Without this (calling each separately) real real	1m0.753s, with:	1m0.478s
func processLine(record []string, header []string, altAlleles []byte, emptyField string,
	fieldDelimiter string, minGq float64, results chan<- string, wg *sync.WaitGroup) {

	defer wg.Done()
	var homs [][]string
	var hets [][]string
	var missing [][]string

	if len(header) > firstSampleIdx {
		homs, hets, missing = makeHetHomozygotes(record, header, altAlleles, minGq)
	}
	// // TODO: process 
	for i, alt := range altAlleles {
		// If no sampels are provided, annotate what we can, skipping hets and homs
		if len(header) > firstSampleIdx {
			if len(homs[i]) == 0 && len(hets[i]) == 0 {
				continue
			}
		}

		var output bytes.Buffer
		output.WriteString(record[chromIdx])
		output.WriteString("\t")
		output.WriteString(record[posIdx])
		output.WriteString("\t")
		output.WriteString(record[typeIdx])
		output.WriteString("\t")
		output.WriteString(record[refIdx])
		output.WriteString("\t")
		output.WriteByte(alt)
		output.WriteString("\t")

		if len(hets[i]) == 0 {
			output.WriteString(emptyField)
		} else {
			output.WriteString(strings.Join(hets[i], fieldDelimiter))
		}

		output.WriteString("\t")

		if len(homs[i]) == 0 {
			output.WriteString(emptyField)
		} else {
			output.WriteString(strings.Join(homs[i], fieldDelimiter))
		}

		output.WriteString("\t")

		if len(missing[i]) == 0 {
			output.WriteString(emptyField)
		} else {
			output.WriteString(strings.Join(missing[i], fieldDelimiter))
		}

		output.WriteString("\n")
		results <- output.String()
	}
}

// Explicitly whitelist call types, because this gets stored, could be malicious
func validType(cType string) bool {
  return cType == "SNP" || cType == "INS" || cType == "DEL" || cType == "MULTIALLELIC" || cType == "DENOVO_SNP" || cType == "DENOVO_INS" || cType == "DENOVO_DEL" || cType == "DENOVO_MULTIALLELIC"
}

func gatherAlt(ref string, alleles string, alt map[string]map[string][]byte) []byte {
	if len(alt[ref][alleles]) > 0 {
		return alt[ref][alleles]
	}

  _, ok := alt[ref];

  if !ok {
    alt[ref] = make(map[string][]byte)
  }

  _, ok = alt[ref][alleles]

  if ok {
    return alt[ref][alleles]
  }

  if !strings.Contains(alleles, ",") {
    // 
    alt[ref][alleles] = []byte{alleles[0]}

    return alt[ref][alleles]
  }

  for _, val := range strings.Split(alleles, ",") {
    if val == ref {
      continue
    }

    alt[ref][alleles] = append(alt[ref][alleles], val[0])
  }

  return alt[ref][alleles]
}

func makeHetHomozygotes(fields []string, header []string, altAlleles []byte, minGQ float64) ([][]string, [][]string, [][]string) {
	missing := make([][]string, len(altAlleles), len(altAlleles))
	het := make([][]string, len(altAlleles), len(altAlleles))
	hom := make([][]string, len(altAlleles), len(altAlleles))

	SAMPLES: for i:= firstSampleIdx; i < len(fields); i+=2 {
		if fields[refIdx] == fields[i] {
			continue
		}

		if fields[i] == "N" {
			appendMissing(len(altAlleles), header[i], missing)
			continue
		}

		conf, err := strconv.ParseFloat(fields[i + 1], 64)

		if err != nil {
			log.Printf("%s:%s: %s genotype invalid confidence %s", fields[chromIdx], fields[posIdx], header[i], fields[i + 1])
			appendMissing(len(altAlleles), header[i], missing)
			continue
		}

		if conf < minGQ {
			appendMissing(len(altAlleles), header[i], missing)
			continue
		}

		//Fast path for homozygotes
		//avoid calculating hash; may be faster for small sets
		//small hash linear search golang
		//GC somewhat depleted in human http://blog.kokocinski.net/index.php/gc-content-of-human-chromosomes?blog=2
		if fields[i][0] == 'A' || fields[i][0] == 'T' || fields[i][0] == 'C' ||
		fields[i][0] == 'G' || fields[i][0] == 'D' || fields[i][0] == 'I' {
			for altIndex, oAlt  := range altAlleles {
				if fields[i][0] == oAlt || (oAlt == '-' && fields[i][0] == 'D') || (oAlt == '+' && fields[i][0] == 'I') {
					hom[altIndex] = append(hom[altIndex], header[i])
					continue SAMPLES
				}
			}

			appendMissing(len(altAlleles), header[i], missing)
			log.Printf("%s:%s: %s genotype %s not in Alleles", fields[chromIdx], fields[posIdx], header[i], fields[i])
			continue
		}

		// Heterozygote
		iupacArr, ok := iupac[fields[i][0]];

		if !ok {
			log.Printf("%s:%s: %s genotype %s not IUPAC", fields[chromIdx], fields[posIdx], header[i], fields[i])
			appendMissing(len(altAlleles), header[i], missing)
			continue
		}

		IUPAC: for _, tAlt := range iupacArr {
			// altAlleles do not include the reference
			// Explicitly we expect all alleles to be 1 char in the input file
			if tAlt == fields[refIdx][0] {
				continue
			}

			// First check that the iupac code makes sense given the present alleles
			// It could be that (in low-enough quality sites) the code doesn't actually correspond
			// to the present alleles
			for altIndex, oAlt  := range altAlleles {
				if tAlt == oAlt {
					het[altIndex] = append(het[altIndex], header[i])
					continue IUPAC
				}
			}

			// If we're here, a sample has a genotype that isn't in Alleles 
			// Since genotype error add it to missingGenos, remove from het list if present
			for altIndex, _ := range altAlleles {
				// If previously added to hets, it would be the last item
				if len(het[altIndex]) > 0 && het[altIndex][len(het[altIndex]) - 1] == header[i]  {
					het[altIndex] = het[altIndex][: len(het[altIndex]) - 1]
				}

				missing[altIndex] = append(missing[altIndex], header[i])
			}

			log.Printf("%s:%s: %s genotype %s not in Alleles", fields[chromIdx], fields[posIdx], header[i], fields[i])
			break IUPAC
		}
	}

	return hom, het, missing
}

func appendMissing(numAlt int, sampleName string, arr [][]string) {
	for i := 0; i < numAlt; i++ {
		arr[i] = append(arr[i], sampleName)
	}
}


func normalizeSampleNames(header []string) {
	re := regexp.MustCompile(`[^a-zA-Z0-9_-]`)

	for i := firstSampleIdx; i < len(header); i+= 2 {
		header[i] = re.ReplaceAllString(header[i], "_")
	}
}
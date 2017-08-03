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
  "strconv"
  "github.com/akotlar/sequtils/parse"
)

const concurrency int = 3
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

func main() {
	config := setup(nil)

	inFh := os.Stdin

	// Close when the functin returns
	defer inFh.Close()

	reader := bufio.NewReader(inFh)

	readSnp(config, reader, func(row string) {fmt.Print(row)})
}

func readSnp(config *Config, reader *bufio.Reader, resultFunc func(row string)) {
	// Read buffer
	workQueue := make(chan string, 100)
	complete := make(chan bool)
	// Write buffer
	results := make(chan string, 100)
	var wg sync.WaitGroup

	endOfLineByte, numChars, headerLine, err := parse.FindEndOfLine(reader, "")

  if err != nil {
    log.Fatal(err)
  }

  header := strings.Split(headerLine[:len(headerLine) - numChars], "\t")

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

			workQueue <- row[:len(row) - numChars];
		}

		// Close the channel so everyone reading from it knows we're done.
		close(workQueue)
	}()

	wg.Add(1)
	go func() {
		defer wg.Done()
		for line := range results {
			resultFunc(line)
		}
	}()

	// Now read them all off, concurrently.
	for i := 0; i < concurrency; i++ {
		go processLine(header, config.emptyField, config.fieldDelimiter, config.minGq, workQueue, results, complete)
	}

	// Wait for everyone to finish.
	for i := 0; i < concurrency; i++ {
		<-complete
	}

	close(results)

	wg.Wait()
}

func processLine(header []string, emptyField string, fieldDelimiter string, minGq float64, queue chan string, results chan string, complete chan bool) {
	alleleCache := make(map[byte]map[string][]string)
	var altAlleles []string

	for line := range queue {
		record := strings.Split(line, "\t")

		if !validType(record[typeIdx]) {
      continue
    }

	  var homs [][]string
	  var hets [][]string
	  var missing [][]string

	  // Hit rate for this is quite high; runs only 12 times for all SNPs
	  // and N times for indels, multiallelics, but those are rare, and do repeat
	  altAlleles = gatherAlt(record[refIdx][0], record[altIdx], alleleCache)

	  if len(header) > firstSampleIdx {
	    homs, hets, missing = makeHetHomozygotes(record, header, altAlleles, minGq)
	  }

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
	    output.WriteString(alt)
	    output.WriteString("\t")

      if(len(altAlleles) > 1) {
        output.WriteRune(parse.NotTrTv)
      } else {
        output.WriteRune(parse.GetTrTv(record[refIdx], alt))
      }

      output.WriteString("\t")

	    if len(hets) == 0 || len(hets[i]) == 0 {
	      output.WriteString(emptyField)
	    } else {
	      output.WriteString(strings.Join(hets[i], fieldDelimiter))
	    }

	    output.WriteString("\t")

	    if len(homs) == 0 || len(homs[i]) == 0 {
	      output.WriteString(emptyField)
	    } else {
	      output.WriteString(strings.Join(homs[i], fieldDelimiter))
	    }

	    output.WriteString("\t")

	    if len(missing) == 0 || len(missing[i]) == 0 {
	      output.WriteString(emptyField)
	    } else {
	      output.WriteString(strings.Join(missing[i], fieldDelimiter))
	    }

	    output.WriteString("\n")
	    results <- output.String()
	  }
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

func makeHetHomozygotes(fields []string, header []string, altAlleles []string, minGQ float64) ([][]string, [][]string, [][]string) {
  missing := make([][]string, len(altAlleles), len(altAlleles))
  het := make([][]string, len(altAlleles), len(altAlleles))
  hom := make([][]string, len(altAlleles), len(altAlleles))

  SAMPLES: for i:= firstSampleIdx; i < len(fields); i+=2 {
    if fields[refIdx] == fields[i] {
      continue
    }

    if fields[i] == "N" {
      parse.AppendMissing(len(altAlleles), header[i], missing)
      continue
    }

    conf, err := strconv.ParseFloat(fields[i + 1], 64)

    if err != nil {
      log.Printf("%s:%s: %s genotype invalid confidence %s", fields[chromIdx], fields[posIdx], header[i], fields[i + 1])
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
      for altIndex, oAlt  := range altAlleles {
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
    iupacArr, ok := iupac[fields[i][0]];

    if !ok {
      log.Printf("%s:%s: %s genotype %s not IUPAC", fields[chromIdx], fields[posIdx], header[i], fields[i])
      parse.AppendMissing(len(altAlleles), header[i], missing)
      continue
    }

    IUPAC: for _, tAlt := range iupacArr {
      if tAlt == fields[refIdx][0] {
        continue
      }

      // First check that the iupac code makes sense given the present alleles
      // It could be that (in low-enough quality sites) the code doesn't actually correspond
      // to the present alleles
      for altIndex, oAlt  := range altAlleles {
        if tAlt == oAlt[0] {
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
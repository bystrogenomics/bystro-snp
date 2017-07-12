package main

import (
	"strings"
	"testing"
	// "fmt"
	"bufio"
	// "log"
	// "io"
	// "sync"
)

func TestKeepFlagsTrue(t *testing.T) {
	args := []string{
		"--inPath", "/path/to/file",
		"--errPath", "/path/to/err",
		"--emptyField", ".",
		"--fieldDelimiter", "&",
	}

	config := setup(args)

  if config.inPath != "/path/to/file" || config.errPath != "/path/to/err" {
  	t.Error("NOT OK: parse inPath and errPath args")
  }

  if config.emptyField != "." || config.fieldDelimiter != "&" {
  	t.Error("NOT OK: parse emptyField and fieldDelimiter args")
  }
}

func TestValidType(t *testing.T) {
  if !validType("SNP") {
    t.Error("Couldn't process SNP")
  }

  if !validType("INS") {
    t.Error("Couldn't process SNP")
  }

  if !validType("DEL") {
    t.Error("Couldn't process SNP")
  }

  if !validType("MULTIALLELIC") {
    t.Error("Couldn't process SNP")
  }

  if !validType("DENOVO_MULTIALLELIC") {
    t.Error("Couldn't process SNP")
  }

  if !validType("DENOVO_SNP") {
    t.Error("Couldn't process SNP")
  }

  if !validType("DENOVO_INS") {
    t.Error("Couldn't process SNP")
  }

  if !validType("DENOVO_DEL") {
    t.Error("Couldn't process SNP")
  }
}

func TestGatherAlt(t *testing.T) {
  altCache := make(map[byte]map[string][]string)

  alleles := gatherAlt('A', "T", altCache)

  if alleles[0] != "T" {
  	t.Error("Expected only alt to be T", alleles)
  }

  if altCache['A']["T"][0] != "T" {
  	t.Error("Expected altCache to contain T alleles for ref:A alleles:T", alleles)
  }

  alleles = gatherAlt('A', "A,T", altCache)

  if alleles[0] != "T" {
  	t.Error("Expected only alt to be T", alleles)
  }

  if altCache['A']["A,T"][0] != "T" {
  	t.Error("Expected altCache to contain T alleles for ref:A alleles:A,T", alleles)
  }

  for _, ref := range []string{"C","T","G"} {
  	if len(altCache[ref[0]]["A,T"]) != 0 {
	  	t.Error("Shouldn't have cached ref", ref, alleles)
  	}
  }

  alleles = gatherAlt('A', "A,T,C", altCache)

  if alleles[0] != "T" || alleles[1] != "C" {
  	t.Error("Expected multiallelic alts to be T first then C", alleles)
  }

  if altCache['A']["A,T,C"][0] != "T" || altCache['A']["A,T,C"][1] != "C" {
  	t.Error("Expected altCache to contain T and C alleles for ref:A alleles:A,T,C", alleles)
  }

  alleles = gatherAlt('A', "A,-9", altCache)

  if alleles[0] != "-9" {
  	t.Error("Expected alt to be - ", alleles)
  }

  if altCache['A']["A,-9"][0] != "-9" {
  	t.Error("Expected altCache to contain - allele for ref:A alleles:A,-9", alleles)
  }

  alleles = gatherAlt('A', "A,+AAT", altCache)

  if alleles[0] != "+AAT" {
  	t.Error("Expected alt to be +", alleles)
  }

  if altCache['A']["A,+AAT"][0] != "+AAT" {
  	t.Error("Expected altCache to contain + allele for ref:A alleles:A,+AAT", alleles)
  }

  alleles = gatherAlt('A', "A,+AATAACCCTTGGGG", altCache)

  if alleles[0] != "+AATAACCCTTGGGG" {
  	t.Error("Expected alt to be +", alleles)
  }

  if altCache['A']["A,+AATAACCCTTGGGG"][0] != "+AATAACCCTTGGGG" {
  	t.Error("Expected altCache to contain + allele for ref:A alleles:A,+AATAACCCTTGGGG", alleles)
  }

  alleles = gatherAlt('A', "A,-9,+AATAACCCTTGGGG", altCache)

  if alleles[0] != "-9" || alleles[1] != "+AATAACCCTTGGGG" {
  	t.Error("Expected alt to be - and + for ref:A, alleles: A,-9,+AATAACCCTTGGGG", alleles)
  }

  if altCache['A']["A,-9,+AATAACCCTTGGGG"][0] != "-9"  || altCache['A']["A,-9,+AATAACCCTTGGGG"][1] != "+AATAACCCTTGGGG" {
  	t.Error("Expected altCache to contain + allele for ref:A alleles:A,-9,+AATAACCCTTGGGG", alleles)
  }

  alleles = gatherAlt('T', "A,-9,+AATAACCCTTGGGG", altCache)

  if alleles[0] != "A" || alleles[1] != "-9" || alleles[2] != "+AATAACCCTTGGGG" {
  	t.Error("Expected alt to be T, -, and + for ref:T, alleles: A,-9,+AATAACCCTTGGGG", alleles)
  }

  if altCache['T']["A,-9,+AATAACCCTTGGGG"][0] != "A" || altCache['T']["A,-9,+AATAACCCTTGGGG"][1] != "-9" || altCache['T']["A,-9,+AATAACCCTTGGGG"][2] != "+AATAACCCTTGGGG" {
  	t.Error("Expected altCache to contain + allele for ref:T alleles:A,-9,+AATAACCCTTGGGG", alleles)
  }

  alleles = gatherAlt('G', "A,C", altCache)

  if alleles[0] != "A" || alleles[1] != "C" {
  	t.Error("Expected alt to be A, C for ref:T, alleles: A,C", alleles)
  }

  if altCache['G']["A,C"][0] != "A" || altCache['G']["A,C"][1] != "C" {
  	t.Error("Expected altCache to contain + allele for ref:T alleles:A,-9,+AATAACCCTTGGGG", alleles)
  }
}

func TestProcessBiAllelicLine(t *testing.T) {
	record := strings.Join([]string{"chr1", "10000", "C", "A,T", "1,100", "SNP", "W", "1", "A", "1", "T", "1"},"\t")
	header := strings.Join([]string{"Fragment", "Position", "Reference", "Alleles", "Allele_counts", "Type", "Sample1", " ", "Sample2", " ", "Sample3", " "}, "\t")

	lines := header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", minGq: .95}

	i := 0
	readSnp(&config, reader, func(row string) {
			record := strings.Split(row[:len(row) - 1], "\t")

		if record[0] != "chr1" || record[1] != "10000" || record[2] != "SNP" || record[3] != "C" {
			t.Error("Expect all rows to have the same chr:pos, reference, and type (SNP) in biallelic SNP")
		}

    if record[5] != "0" {
      t.Error("Expect bi-allelic snp to not be tr or tv with value 0")
    }

		if i == 0 {
      if record[4] != "A" {
        t.Error("Expected allele 1 to be A", record)
      }

			if record[6] != "Sample1" {
				t.Error("Expect Sample1 to be shown as het for allele A", record)
			}

			if record[7] != "Sample2" {
				t.Error("Expect Sample2 to be shown as homozygous for allele A", record)
			}

			if record[8] != "!" {
				t.Error("Expect no missing genotypes", record)
			}
		}

		if i == 1 {
      if record[4] != "T" {
        t.Error("Expected allele 2 to be T", record)
      }

			if record[6] != "Sample1" {
				t.Error("Expect Sample1 to be shown as het for allele T", record)
			}

			if record[7] != "Sample3" {
				t.Error("Expect Sample3 to be shown as het for allele T", record)
			}

			if record[8] != "!" {
				t.Error("Expect no missing genotypes", record)
			}
		}

		i++
	})

	if i != 2 {
		t.Error("Expect 2 alleles to be parsed from biallelic SNP, got: ", i)
	}
}

func TestProcessBiAllelicLineWithGenotypingError(t *testing.T) {
	record := strings.Join([]string{"chr1", "10000", "C", "A,T", "1,100", "SNP", "K", "1", "A", "1", "T", "1"}, "\t")
	header := strings.Join([]string{"Fragment", "Position", "Reference", "Alleles", "Allele_counts", "Type", "Sample1", " ", "Sample2", " ", "Sample3", " "}, "\t")

	lines := header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", minGq: .95}

	i := 0
	readSnp(&config, reader, func(row string) {
		record := strings.Split(row[:len(row) - 1], "\t")

		if record[0] != "chr1" || record[1] != "10000" || record[2] != "SNP" || record[3] != "C" {
			t.Error("Expect all rows to have the same chr:pos, reference, and type (SNP) in biallelic SNP")
		}

    if record[5] != "0" {
      t.Error("Expect bi-allelic snp to not be tr or tv with value 0")
    }

		if i == 0 {
      if record[4] != "A" {
        t.Error("Expected allele 1 to be A", record)
      }

			if record[6] != "!" {
				t.Error("Expect Sample1 to be shown as het for allele A", record)
			}

			if record[7] != "Sample2" {
				t.Error("Expect Sample2 to be shown as homozygous for allele A", record)
			}

			if record[8] != "Sample1" {
				t.Error("Expect Sample1 to be missing for first allele, because miscalled relative to Alleles", record)
			}
		}

		if i == 1 {
      if record[4] != "T" {
        t.Error("Expected allele 2 to be T", record)
      }

			if record[6] != "!" {
				t.Error("Expect Sample1 to be shown as het for allele T", record)
			}

			if record[7] != "Sample3" {
				t.Error("Expect Sample3 to be shown as het for allele T", record)
			}

			if record[8] != "Sample1" {
				t.Error("Expect Sample1 to be missing for first allele, because miscalled relative to Alleles", record)
			}
		}

		i++
	})

	if i != 2 {
		t.Error("Expect 2 alleles to be parsed from biallelic SNP")
	}
}

func TestProcessBiAllelicLineWithLowCoverageError(t *testing.T) {
	record := strings.Join([]string{"chr1", "10000", "C", "A,T", "1,100", "SNP", "K", "1", "A", "1", "T", ".9"},"\t")
	header := strings.Join([]string{"Fragment", "Position", "Reference", "Alleles", "Allele_counts", "Type", "SampleA", " ", "Sample2", " ", "SampleB", " "}, "\t")

	lines := header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", minGq: .95}

	i := 0
	readSnp(&config, reader, func(row string) {
		record := strings.Split(row[:len(row) - 1], "\t")

		if record[0] != "chr1" || record[1] != "10000" || record[2] != "SNP" || record[3] != "C" {
			t.Error("Expect all rows to have the same chr:pos, reference, and type (SNP) in biallelic SNP")
		}

    if record[5] != "0" {
      t.Error("Expect bi-allelic snp to not be tr or tv with value 0")
    }

		if i == 0 {
      if record[4] != "A" {
        t.Error("Expected allele 1 to be A", record)
      }

			if record[6] != "!" {
				t.Error("Expect Sample1 to be shown as het for allele A", record)
			}

			if record[7] != "Sample2" {
				t.Error("Expect Sample2 to be shown as homozygous for allele A", record)
			}

			if record[8] != "SampleA;SampleB" {
				t.Error(`Expect Sample1 and Sample3 to be missing for first allele, because
					Samples 1 miscalled relative to Alleles, and Sample3 lower than requested .95 confidence`, record)
			}
		}

		// won't be reached; because allele T has no valid alleles
		// if i == 1 {
		// }

		i++
	})

	if i != 1 {
		t.Error("Expect 2 alleles to be parsed from biallelic SNP")
	}
}

func TestProcessMultiallelicLine(t *testing.T) {
	header := strings.Join([]string{"Fragment", "Position", "Reference", "Alleles", "Allele_counts", "Type",
	"Sample1", " ", "Sample2", " ", "Sample3", " ", "Sample4", " ", "Sample5", " ", "Sample6", " ", "Sample7", " ", "Sample8", " ",
	"Sample9", " ", "Sample10", " ", "Sample11", " "}, "\t")
	record := strings.Join([]string{"chr1", "10000", "C", "C,T,+AATC,-9", "1,100", "MULTIALLELIC",
	"E", "1", "H", "1", "D", "1", "I", "1", "C", "1", "Y", "1", "R", "1", "K", ".9", "T", ".94", "T", ".95", "Y", ".949"}, "\t")
	//1        2         3         4         5         6         7         8          9           10          11

	lines := header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", minGq: .95}

	i := 0
	readSnp(&config, reader, func(row string) {
		record := strings.Split(row[:len(row) - 1], "\t")

		if record[0] != "chr1" || record[1] != "10000" || record[2] != "MULTIALLELIC" || record[3] != "C" {
			t.Error("Expect all rows to have the same chr:pos, reference, and type (MULTIALLELIC) in MULTIALLELIC")
		}

    if record[5] != "0" {
      t.Error("Expect multiallelics to be set to 0 trTv", record)
    }

		// Allele 1 is reference, will not be put into a row
		// Likewise, Sample 5 is reference homozygote, so will not appear in any row
		// the T alelles
		if i == 0 {
      if record[4] != "T" {
        t.Error("Expect allele 1 to be C", record)
      }
			// Sample6 has genotype Y == [C,T], so that is a valid het
			// Sample7 has genotype R == [A,G], so that is invalid (missing)
			// Sample8 has genotype K == [G,T], so that is invalid (missing)
			// Sample9 has genotype T, which is valid homozygote, but it has too low confidence, so missing in all rows
			// Sample10 has genotype T, which is valid homozygote, and has just high enough confidence
			// Sample11 has genotype Y, which is valid het, but confidence is too low
			if record[6] != "Sample6" {
				t.Error("Expect Sample6 to be shown as het for allele T", record)
			}

			if record[7] != "Sample10" {
				t.Error("Expect Sample10 to be shown as homozygous for allele T", record)
			}

			if record[8] != "Sample7;Sample8;Sample9;Sample11" {
				t.Error("Expect Sample7;Sample8;Sample9;Sample11 missing on all lines", record)
			}
		}

		// The +AATC
		if i == 1 {
      if record[4] != "+AATC" {
        t.Error("Expect allele 1 to be +AATC", record)
      }

			if record[6] != "Sample2" {
				t.Error("Expect Sample2 to be shown as het for allele +AATC", record)
			}

			if record[7] != "Sample4" {
				t.Error("Expect Sample3 to be shown as homozygous for allele +AATC", record)
			}

			if record[8] != "Sample7;Sample8;Sample9;Sample11" {
				t.Error("Expect Sample7;Sample8;Sample9;Sample11 missing on all lines", record)
			}
		}

		// The -9
		if i == 2 {
      if record[4] != "-9" {
        t.Error("Expect allele 1 to be -9", record)
      }

			if record[6] != "Sample1" {
				t.Error("Expect Sample1 to be shown as het for allele -9", record)
			}

			if record[7] != "Sample3" {
				t.Error("Expect Sample3 to be shown as homozygous for allele -9", record)
			}

			if record[8] != "Sample7;Sample8;Sample9;Sample11" {
				t.Error("Expect Sample7;Sample8;Sample9;Sample11 missing on all lines", record)
			}
		}

		i++
	})

	if i != 3 {
		t.Error("Expect 3 alleles to be parsed from this multiallelic", header, record)
	}
}

func TestProcessSingleLine(t *testing.T) {
	header := strings.Join([]string{"Fragment", "Position", "Reference", "Alleles", "Allele_counts", "Type",
	"Sample1", " ", "Sample2", " ", "Sample3", " ", "Sample4", " ", "Sample5", " ", "Sample6", " ", "Sample7", " ", "Sample8", " ",
	"Sample9", " ", "Sample10", " ", "Sample11", " ", "Sample12", " "}, "\t")
	record := strings.Join([]string{"chr11", "12000", "C", "C,T", "1,100", "SNP",
	"C", "1", "T", "1", "Y", "1", "Y", "1", "C", "1", "C", "1", "C", "1", "K", "1", "Y", ".9", "R", ".94", "G", ".95", "W", ".949"}, "\t")
	//1        2         3         4         5         6         7         8          9           10          11

	lines := header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", minGq: .95}

	i := 0
	readSnp(&config, reader, func(row string) {
		record := strings.Split(row[:len(row) - 1], "\t")

		if record[0] != "chr11" || record[1] != "12000" || record[2] != "SNP" || record[3] != "C" {
			t.Error("Expect all rows to have the same chr:pos, reference, and type (SNP) in biallelic SNP", record)
		}

    if record[4] != "T" {
        t.Error("Expected allele 1 to be T", record)
    }

    if record[5] != "1" {
      t.Error("Expect SNP C->T to be a transition with value 1")
    }

		// Sample6 has genotype Y == [C,T], so that is a valid het
		// Sample7 has genotype R == [A,G], so that is invalid (missing)
		// Sample8 has genotype K == [G,T], so that is invalid (missing)
		// Sample9 has genotype T, which is valid homozygote, but it has too low confidence, so missing in all rows
		// Sample10 has genotype T, which is valid homozygote, and has just high enough confidence
		// Sample11 has genotype Y, which is valid het, but confidence is too low
		if record[6] != "Sample3;Sample4" {
			t.Error("Expect Sample6 to be shown as het for allele T", record)
		}

		if record[7] != "Sample2" {
			t.Error("Expect Sample2 to be shown as homozygous for allele for simple SNP", record)
		}

		if record[8] != "Sample8;Sample9;Sample10;Sample11;Sample12" {
			t.Error("Expect Sample8;Sample9;Sample10;Sample11;Sample12 missing on all lines", record)
		}

		i++
	})

	if i != 1 {
		t.Error("Expect 3 alleles to be parsed from this multiallelic", header, record)
	}
}

func TestProcessMultiallelicSnpLine(t *testing.T) {
	header := strings.Join([]string{"Fragment", "Position", "Reference", "Alleles", "Allele_counts", "Type",
	"Sample1", " ", "Sample2", " ", "Sample3", " ", "Sample4", " ", "Sample5", " ", "Sample6", " ", "Sample7", " ", "Sample8", " ",
	"Sample9", " ", "Sample10", " ", "Sample11", " ", "Sample12", " ", "Sample13", " ", "Sample14", " "}, "\t")
	record := strings.Join([]string{"chr1", "10000", "A", "C,T,A,G", "1,100", "MULTIALLELIC",
	"A", "1", "C", "1", "T", "1", "G", "1", "R", "1", "Y", "1", "G", ".5", "S", "1", "T", "1", "W", "1", "K", "1", "M", "1", "C", ".9", "E", ".96"}, "\t")
	//1        2         3         4         5         6         7         8         9         10        11        12        13          14

	lines := header + "\n" + record	+ "\n"

	reader := bufio.NewReader(strings.NewReader(lines))

	config := Config{emptyField: "!", fieldDelimiter: ";", minGq: .95}

	i := 0
	readSnp(&config, reader, func(row string) {
		record := strings.Split(row[:len(row) - 1], "\t")

		if record[0] != "chr1" || record[1] != "10000" || record[2] != "MULTIALLELIC" || record[3] != "A"{
			t.Error("Expect all rows to have the same chr:pos, reference, and type (MULTIALLELIC) in MULTIALLELIC", record)
		}

    if record[5] != "0" {
      t.Error("Expect multiallelic to have trTv 0")
    }
		// Allele 3 is reference, will not be put into a row
		if i == 0 {
			// A -> C
      if record[4] != "C" {
        t.Error("Expect first allele to be C")
      }
			// Sample2 has genotype C, homozygote
			// Sample6 has genotype Y == [C,T], het (will also be het for allele 2, T)
			// Sample8 has genotype S == [G,C], het (will also be het for allele 4, G)
			// Sample12 has genotype M, het
			// Sample7 has G, which doesn't match, but it is missing because .5 < .95 confidence
			// Sample13 has C, but .9 < .95 confidence, so missing in all alleles
			// Sample14 has an E, which doesn't exist in the allele list
			if record[6] != "Sample6;Sample8;Sample12" {
				t.Error("Expect Sample6 to be shown as het for allele T", record)
			}

			if record[7] != "Sample2" {
				t.Error("Expect Sample10 to be shown as homozygous for allele T", record)
			}

			if record[8] != "Sample7;Sample13;Sample14" {
				t.Error("Expect Sample7;Sample8;Sample9;Sample11 missing on all lines", record)
			}
		}

		if i == 1 {
			//A -> T
      if record[4] != "T" {
        t.Error("Expect second allele to be T")
      }
			//Sample3 is T, homozygote
			//Sample6 is Y == [C, T], so it is heterozygous for T (also for C, allele 1)
			//Sample9 is also T, homozygous  
			//Sample10 is W == [A, T], so is het
			//Sample11 is K == [G, T], so is het (also het for allele 4, G)
			if record[6] != "Sample6;Sample10;Sample11" {
				t.Error("Expect Sample6 to be shown as het for allele T", record)
			}

			if record[7] != "Sample3;Sample9" {
				t.Error("Expect Sample10 to be shown as homozygous for allele T", record)
			}

			if record[8] != "Sample7;Sample13;Sample14" {
				t.Error("Expect Sample7;Sample8;Sample9;Sample11 missing on all lines", record)
			}
		}

		i++
	})

	if i != 3 {
		t.Error("Expect 3 alleles to be parsed from this multiallelic", header, record)
	}
}

// func TestFindEndOfLineChar (t *testing.T) {
// 	header := "Testing 1 2 3\r"
// 	reader := bufio.NewReader
// }
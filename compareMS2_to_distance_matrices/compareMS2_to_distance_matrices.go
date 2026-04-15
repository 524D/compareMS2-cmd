package main

// SPDX-License-Identifier: MIT
// Copyright 2022 Magnus Palmblad
// Copyright 2022 Rob Marissen.
// Go conversion with help from GitHub Copilot with Claude Sonnet 4.5

// compareMS2_to_distance_matrices - collects output from compareMS2 to distance matrices
// Contact information: n.m.palmblad@lumc.nl

import (
	"bufio"
	"encoding/json"
	"flag"
	"fmt"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"unicode"
)

// SampleToSpecies represents a species with its associated samples
type SampleToSpecies struct {
	SpeciesName string
	SampleNames []string
	SpeciesUsed bool
}

// Species contains all species data
type Species struct {
	Data []SampleToSpecies
}

// DistanceData holds distance values and sample counts
type DistanceData struct {
	Value   float64
	Samples int
}

// QCData holds QC values and sample counts
type QCData struct {
	Value   float64
	Samples int
}

// ComparisonResult for JSON parsing
type ComparisonResult struct {
	DatasetA    string  `json:"datasetA"`
	DatasetB    string  `json:"datasetB"`
	SetDistance float64 `json:"setDistance"`
	DatasetAQC  float64 `json:"datasetAQC"`
	DatasetBQC  float64 `json:"datasetBQC"`
}

// DistanceMatrixOutput for JSON output
type DistanceMatrixOutput struct {
	DistanceMatrix struct {
		Title     string      `json:"title"`
		Cutoff    float64     `json:"cutoff"`
		Species   []SpeciesQC `json:"species"`
		Distances [][]float64 `json:"distances"`
	} `json:"distanceMatrix"`
}

// SpeciesQC for JSON output
type SpeciesQC struct {
	Name string  `json:"name"`
	QC   float64 `json:"qc"`
}

// Add a species to the species list
func (s *Species) addSpecies(speciesName string) {
	s.Data = append(s.Data, SampleToSpecies{
		SpeciesName: speciesName,
		SampleNames: make([]string, 0),
		SpeciesUsed: false,
	})
}

// Find species index by name, add if not found
func (s *Species) speciesNameToIdx(speciesName string) int {
	for i := range s.Data {
		if s.Data[i].SpeciesName == speciesName {
			return i
		}
	}
	s.addSpecies(speciesName)
	return len(s.Data) - 1
}

// Add a sample to a species
func (s *Species) addSampleSpecies(sampleName, speciesName string) int {
	idx := s.speciesNameToIdx(speciesName)
	s.Data[idx].SampleNames = append(s.Data[idx].SampleNames, sampleName)
	return idx
}

// Get basename from path (works for both Unix and Windows paths)
func basename(path string) string {
	// Use filepath.Base for cross-platform support
	base := filepath.Base(path)
	return base
}

// Remove extension from filename
func removeExtension(filename string) string {
	ext := filepath.Ext(filename)
	if ext != "" {
		return filename[:len(filename)-len(ext)]
	}
	return filename
}

// Map sample name to species index
func (s *Species) sampleNameToSpeciesIndex(sampleName string) int {
	base := basename(sampleName)

	// Search for existing species with this sample
	for i := range s.Data {
		for _, sample := range s.Data[i].SampleNames {
			if sample == base {
				s.Data[i].SpeciesUsed = true
				return i
			}
		}
	}

	// Not found, create new species from sample name
	speciesName := removeExtension(base)
	idx := len(s.Data)
	s.Data = append(s.Data, SampleToSpecies{
		SpeciesName: speciesName,
		SampleNames: []string{base},
		SpeciesUsed: true,
	})
	return idx
}

// Read sample-to-species mapping file
func readSampleToSpecies(filename string, species *Species) error {
	file, err := os.Open(filename)
	if err != nil {
		return fmt.Errorf("error opening sample-species mapping file %s: %v", filename, err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}

		parts := strings.Split(line, "\t")
		if len(parts) >= 2 {
			sampleName := strings.TrimSpace(parts[0])
			speciesName := strings.TrimSpace(parts[1])
			species.addSampleSpecies(sampleName, speciesName)
		}
	}

	return scanner.Err()
}

// Read file list
func readFileList(filename string) ([]string, error) {
	file, err := os.Open(filename)
	if err != nil {
		return nil, fmt.Errorf("error opening filename list %s: %v", filename, err)
	}
	defer file.Close()

	var files []string
	scanner := bufio.NewScanner(file)
	for scanner.Scan() {
		line := strings.TrimSpace(scanner.Text())
		if line != "" {
			files = append(files, line)
		}
	}

	if err := scanner.Err(); err != nil {
		return nil, err
	}

	return files, nil
}

// Get number of species with data
func (s *Species) getNumSpecies() int {
	count := 0
	for i := range s.Data {
		if s.Data[i].SpeciesUsed {
			count++
		}
	}
	return count
}

// Calculate distance matrix index for lower triangular matrix
// where x < y
func distanceIndex(x, y int) int {
	return ((y-1)*y)/2 + x
}

// Create NEXUS format output
func createNexus(outputStem string, distances map[int]DistanceData, species *Species, metric int) error {
	filename := outputStem + "_distance_matrix.nexus"
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", filename, err)
	}
	defer file.Close()

	fmt.Println("writing distance matrix in NEXUS format...")

	nSpecies := species.getNumSpecies()

	fmt.Fprintln(file, "#NEXUS")
	fmt.Fprintln(file, "[                                                                              ]")
	fmt.Fprintln(file, "[ NEXUS file generated by compareMS2_to_distance_matrices from compareMS2 data ]")
	fmt.Fprintln(file, "[ compareMS2 and compareMS2_to_distance_matrices are (c) Magnus Palmblad 2010- ]")
	fmt.Fprintln(file, "[                                                                              ]")
	fmt.Fprintln(file, "[ note that compareMS2_to_distance_matrices outputs a full matrix and that the ]")
	fmt.Fprintln(file, "[ upper (right) and lower (left) triangular matrices can be slightly different ]")
	fmt.Fprintln(file, "[                                                                              ]")
	fmt.Fprintf(file, "[ distances are based on metric %d                                              ]\n", metric)
	fmt.Fprintln(file, "[                                                                              ]")
	fmt.Fprintln(file, "BEGIN taxa;")
	fmt.Fprintf(file, "   DIMENSIONS ntax=%d;\n", nSpecies)
	fmt.Fprint(file, "   TAXLABELS")

	for i := range species.Data {
		if species.Data[i].SpeciesUsed {
			fmt.Fprintf(file, "   %s", species.Data[i].SpeciesName)
		}
	}
	fmt.Fprintln(file, ";")
	fmt.Fprintln(file, "END;")
	fmt.Fprintln(file, "")
	fmt.Fprintln(file, "BEGIN distances;")
	fmt.Fprintf(file, "   DIMENSIONS ntax=%d;\n", nSpecies)
	fmt.Fprintln(file, "   FORMAT")
	fmt.Fprintln(file, "       triangle=both")
	fmt.Fprintln(file, "       labels=left")
	fmt.Fprintln(file, "     diagonal")
	fmt.Fprintln(file, "     missing=?")
	fmt.Fprintln(file, "   ;")
	fmt.Fprintln(file, "MATRIX")

	for i := range species.Data {
		if species.Data[i].SpeciesUsed {
			anyOut := false
			for j := range species.Data {
				if species.Data[j].SpeciesUsed {
					if !anyOut {
						fmt.Fprintf(file, "%s\t", species.Data[i].SpeciesName)
						anyOut = true
					}

					if i == j {
						fmt.Fprintf(file, "%1.5f\t", 0.0)
					} else if i < j {
						dist := distances[distanceIndex(i, j)].Value
						fmt.Fprintf(file, "%1.5f\t", dist)
					} else {
						// Upper triangular - swap coordinates
						dist := distances[distanceIndex(j, i)].Value
						fmt.Fprintf(file, "%1.5f\t", dist)
					}
				}
			}
			if anyOut {
				fmt.Fprintln(file, "")
			}
		}
	}

	fmt.Fprintln(file, ";")
	fmt.Fprintln(file, "END;")

	return nil
}

// Create MEGA format (version < 12)
func createMega(outputStem string, distances map[int]DistanceData, species *Species, cutoff float64, qcData map[int]QCData) error {
	filename := outputStem + "_distance_matrix.meg"
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", filename, err)
	}
	defer file.Close()

	fmt.Println("writing distance matrix in MEGA version < 12 format...")

	fmt.Fprintln(file, "#mega")
	fmt.Fprintf(file, "TITLE: %s (lower-left triangular matrix, cutoff=%.4f)\n\n", filename, cutoff)

	// Output QC values
	for i := range species.Data {
		if species.Data[i].SpeciesUsed {
			qc := qcData[i]
			avgQC := 0.0
			if qc.Samples > 0 {
				avgQC = qc.Value / float64(qc.Samples)
			}
			fmt.Fprintf(file, "QC\t%s\t%.3f\n", species.Data[i].SpeciesName, avgQC)
		}
	}
	fmt.Fprintln(file, "")

	// Output species names
	for i := range species.Data {
		if species.Data[i].SpeciesUsed {
			fmt.Fprintf(file, "#%s\n", species.Data[i].SpeciesName)
		}
	}
	fmt.Fprintln(file, "\n")

	// Output distance matrix (lower triangular)
	for y := 1; y < len(species.Data); y++ {
		if species.Data[y].SpeciesUsed {
			anyOut := false
			for x := 0; x < y; x++ {
				if species.Data[x].SpeciesUsed {
					dist := distances[distanceIndex(x, y)].Value
					fmt.Fprintf(file, "%1.5f\t", dist)
					anyOut = true
				}
			}
			if anyOut {
				fmt.Fprintln(file, "")
			}
		}
	}

	return nil
}

// Fix taxonomy name for MEGA format
func fixTaxonomyName(name string) string {
	runes := []rune(name)
	if len(runes) == 0 {
		return name
	}

	// First character
	if !unicode.IsLetter(runes[0]) && !unicode.IsDigit(runes[0]) &&
		runes[0] != '-' && runes[0] != '+' && runes[0] != '.' {
		runes[0] = '_'
	}

	// Remaining characters
	for i := 1; i < len(runes); i++ {
		if !unicode.IsLetter(runes[i]) && !unicode.IsDigit(runes[i]) &&
			runes[i] != '-' && runes[i] != '+' && runes[i] != '.' &&
			runes[i] != '_' && runes[i] != '*' && runes[i] != ':' &&
			runes[i] != '(' && runes[i] != ')' && runes[i] != '|' &&
			runes[i] != '\\' && runes[i] != '/' {
			runes[i] = '_'
		}
	}

	return string(runes)
}

// Create MEGA format version 12
func createMega12(outputStem string, distances map[int]DistanceData, species *Species, cutoff float64) error {
	filename := outputStem + "_distance_matrix.meg"
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", filename, err)
	}
	defer file.Close()

	fmt.Println("writing distance matrix in MEGA 12 format...")

	nSpecies := species.getNumSpecies()

	fmt.Fprintln(file, "#mega")
	fmt.Fprintf(file, "!Title: %s;\n", outputStem)
	fmt.Fprintf(file, "!Format DataType=Distance DataFormat=LowerLeft NTaxa=%d;\n", nSpecies)
	fmt.Fprintln(file, "!Description")
	fmt.Fprintf(file, "  No. of Groups : %d\n", nSpecies)
	fmt.Fprintf(file, "  Cutoff : %.4f\n", cutoff)
	fmt.Fprintln(file, ";")
	fmt.Fprintln(file, "")

	// Output numbered species list
	speciesIdx := 1
	for i := range species.Data {
		if species.Data[i].SpeciesUsed {
			fixedName := fixTaxonomyName(species.Data[i].SpeciesName)
			fmt.Fprintf(file, "[%d] #%s\n", speciesIdx, fixedName)
			speciesIdx++
		}
	}
	fmt.Fprintln(file, "")

	// Output column header
	fmt.Fprint(file, "[")
	speciesIdx = 1
	for i := range species.Data {
		if species.Data[i].SpeciesUsed {
			fmt.Fprintf(file, "%9d ", speciesIdx)
			speciesIdx++
		}
	}
	fmt.Fprintln(file, " ]")

	// Output distance matrix
	speciesIdx = 1
	for y := range species.Data {
		if species.Data[y].SpeciesUsed {
			fmt.Fprintf(file, "[%d] ", speciesIdx)
			for x := range species.Data {
				if species.Data[x].SpeciesUsed {
					if x < y {
						// Lower triangular part
						dist := distances[distanceIndex(x, y)].Value
						fmt.Fprintf(file, "%9.5f ", dist)
					} else if x == y {
						// Diagonal - empty
						fmt.Fprint(file, "          ")
					}
					// Upper triangular - nothing
				}
			}
			fmt.Fprintln(file, "")
			speciesIdx++
		}
	}
	fmt.Fprintln(file, "")

	return nil
}

// Create JSON format output
func createJSON(outputStem string, distances map[int]DistanceData, species *Species, cutoff float64, qcData map[int]QCData) error {
	filename := outputStem + "_distance_matrix.json"
	file, err := os.Create(filename)
	if err != nil {
		return fmt.Errorf("error opening output file %s: %v", filename, err)
	}
	defer file.Close()

	fmt.Println("writing distance matrix in JSON format...")

	var output DistanceMatrixOutput
	output.DistanceMatrix.Title = outputStem
	output.DistanceMatrix.Cutoff = cutoff

	// Build species list with QC values
	output.DistanceMatrix.Species = make([]SpeciesQC, 0)
	for i := range species.Data {
		if species.Data[i].SpeciesUsed {
			qc := qcData[i]
			avgQC := 0.0
			if qc.Samples > 0 {
				avgQC = qc.Value / float64(qc.Samples)
			}
			output.DistanceMatrix.Species = append(output.DistanceMatrix.Species, SpeciesQC{
				Name: species.Data[i].SpeciesName,
				QC:   avgQC,
			})
		}
	}

	// Build distance matrix (lower triangular)
	output.DistanceMatrix.Distances = make([][]float64, 0)
	for y := 1; y < len(species.Data); y++ {
		if species.Data[y].SpeciesUsed {
			row := make([]float64, 0)
			for x := 0; x < y; x++ {
				if species.Data[x].SpeciesUsed {
					dist := distances[distanceIndex(x, y)].Value
					row = append(row, dist)
				}
			}
			if len(row) > 0 {
				output.DistanceMatrix.Distances = append(output.DistanceMatrix.Distances, row)
			}
		}
	}

	encoder := json.NewEncoder(file)
	encoder.SetIndent("", "  ")
	if err := encoder.Encode(output); err != nil {
		return fmt.Errorf("error encoding JSON: %v", err)
	}

	return nil
}

// Parse JSON comparison file
func parseComparisonJSON(filename string, compIdx int, species *Species, distances map[int]DistanceData, qcData map[int]QCData) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	fmt.Printf("\nreading pairwise comparison %d (JSON)", compIdx)

	// Try to parse as proper JSON using the ComparisonResult struct
	var result ComparisonResult
	decoder := json.NewDecoder(file)
	if err := decoder.Decode(&result); err != nil {
		return fmt.Errorf("error decoding JSON file %s: %v", filename, err)
	}

	a := species.sampleNameToSpeciesIndex(result.DatasetA)
	b := species.sampleNameToSpeciesIndex(result.DatasetB)

	// Ensure lower triangular: x < y
	x, y := a, b
	if a > b {
		x, y = b, a
	}

	idx := distanceIndex(x, y)
	d := distances[idx]
	d.Value += result.SetDistance
	d.Samples++
	distances[idx] = d

	qcA := qcData[a]
	qcA.Value += result.DatasetAQC
	qcA.Samples++
	qcData[a] = qcA

	qcB := qcData[b]
	qcB.Value += result.DatasetBQC
	qcB.Samples++
	qcData[b] = qcB

	return nil
}

// Parse text comparison file
func parseComparisonTxt(filename string, compIdx int, species *Species, distances map[int]DistanceData, qcData map[int]QCData) error {
	file, err := os.Open(filename)
	if err != nil {
		return err
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	var a, b, x, y int

	for scanner.Scan() {
		line := scanner.Text()
		if line == "" {
			continue
		}

		parts := strings.Split(line, "\t")
		if len(parts) < 2 {
			continue
		}

		field := parts[0]
		value := parts[1]

		switch field {
		case "dataset_A":
			a = species.sampleNameToSpeciesIndex(value)
			fmt.Printf("\nread pairwise comparison %d (%s and ", compIdx, value)

		case "dataset_B":
			b = species.sampleNameToSpeciesIndex(value)
			// Ensure lower triangular
			if a < b {
				x, y = a, b
			} else {
				x, y = b, a
			}
			fmt.Printf("%s)", value)

		case "set_distance":
			dist, _ := strconv.ParseFloat(value, 64)
			idx := distanceIndex(x, y)
			d := distances[idx]
			d.Value += dist
			d.Samples++
			distances[idx] = d

		case "dataset_A_QC":
			qc, _ := strconv.ParseFloat(value, 64)
			qcA := qcData[a]
			qcA.Value += qc
			qcA.Samples++
			qcData[a] = qcA

		case "dataset_B_QC":
			qc, _ := strconv.ParseFloat(value, 64)
			qcB := qcData[b]
			qcB.Value += qc
			qcB.Samples++
			qcData[b] = qcB
		}
	}

	return scanner.Err()
}

func usage() {
	fmt.Println("usage: compareMS2_to_distance_matrices -i <input file> -o <output stem> [options]")
	fmt.Println()
	fmt.Println("Required arguments:")
	fmt.Println("  -i <file>    List of compareMS2 results files to process")
	fmt.Println("  -o <stem>    Output file stem (output filename without extension)")
	fmt.Println()
	fmt.Println("Optional arguments:")
	fmt.Println("  -x <file>    Sample to species mapping file")
	fmt.Println("  -c <value>   Score cutoff (default: 0.80)")
	fmt.Println("  -n           Output in NEXUS format (default if no format specified)")
	fmt.Println("  -m           Output in MEGA format (version < 12)")
	fmt.Println("  -m2          Output in MEGA format (version 12+)")
	fmt.Println("  -J           Output in JSON format")
	fmt.Println("  -h           Display this help message")
	fmt.Println()
	fmt.Println("Note: Multiple format flags can be combined to generate multiple outputs.")
	fmt.Println()
	fmt.Println("Examples:")
	fmt.Println("  compareMS2_to_distance_matrices -i filelist.txt -o results -m2")
	fmt.Println("  compareMS2_to_distance_matrices -i filelist.txt -o results -m2 -J")
}

func main() {
	// Define flags
	inputFile := flag.String("i", "", "List of compareMS2 results files to process")
	outputStem := flag.String("o", "", "Output file stem (output filename without extension)")
	mappingFile := flag.String("x", "", "Sample to species mapping file")
	cutoff := flag.Float64("c", 0.80, "Score cutoff")
	formatNexus := flag.Bool("n", false, "Output in NEXUS format")
	formatMega := flag.Bool("m", false, "Output in MEGA format (version < 12)")
	formatMega12 := flag.Bool("m2", false, "Output in MEGA format (version 12+)")
	formatJSON := flag.Bool("J", false, "Output in JSON format")
	help := flag.Bool("h", false, "Display help message")

	flag.Parse()

	if *help {
		fmt.Println("compareMS2_to_distance_matrices - (c) Magnus Palmblad")
		fmt.Println()
		usage()
		os.Exit(0)
	}

	if *inputFile == "" || *outputStem == "" {
		usage()
		os.Exit(1)
	}

	// Determine output formats
	hasFormat := *formatNexus || *formatMega || *formatMega12 || *formatJSON
	if !hasFormat {
		// Default to NEXUS if no format specified
		*formatNexus = true
	}

	fmt.Print("reading list of compareMS2 results files...")

	// Read comparison file list
	comparisons, err := readFileList(*inputFile)
	if err != nil {
		fmt.Fprintf(os.Stderr, "\n%v\n", err)
		os.Exit(1)
	}

	species := &Species{Data: make([]SampleToSpecies, 0)}

	// Read sample-to-species mapping if provided
	if *mappingFile != "" {
		fmt.Print("done\nreading sample to species file...")
		if err := readSampleToSpecies(*mappingFile, species); err != nil {
			fmt.Fprintf(os.Stderr, "\n%v\n", err)
			os.Exit(1)
		}
	}

	fmt.Println("done")

	// Initialize data structures using maps
	distances := make(map[int]DistanceData)
	qcData := make(map[int]QCData)

	// Process comparison files
	for i, compFile := range comparisons {
		isJSON := strings.HasSuffix(compFile, ".json")

		var err error
		if isJSON {
			err = parseComparisonJSON(compFile, i+1, species, distances, qcData)
		} else {
			err = parseComparisonTxt(compFile, i+1, species, distances, qcData)
		}

		if err != nil {
			fmt.Fprintf(os.Stderr, "\nerror opening MS2compare results file %s: %v\n", compFile, err)
			os.Exit(1)
		}
	}

	fmt.Printf("\ndone (read %d pairwise comparisons)\n", len(comparisons))

	// Compute averages for distances with multiple samples
	// FIXME: different from C - the C code uses n_comparisons as the loop limit,
	// but should iterate over all distance entries. This Go version correctly
	// iterates over the actual distance data.
	for idx, d := range distances {
		if d.Samples > 0 {
			d.Value /= float64(d.Samples)
			distances[idx] = d
		}
	}

	// Generate output in requested formats
	metric := 0 // Not used in current implementation

	if *formatNexus {
		if err := createNexus(*outputStem, distances, species, metric); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			os.Exit(1)
		}
	}

	if *formatMega {
		if err := createMega(*outputStem, distances, species, *cutoff, qcData); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			os.Exit(1)
		}
	}

	if *formatMega12 {
		// When both MEGA and MEGA12 are requested, add suffix to avoid overwriting
		stem := *outputStem
		if *formatMega {
			stem = *outputStem + "_MEGA_12_"
		}
		if err := createMega12(stem, distances, species, *cutoff); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			os.Exit(1)
		}
	}

	if *formatJSON {
		if err := createJSON(*outputStem, distances, species, *cutoff, qcData); err != nil {
			fmt.Fprintf(os.Stderr, "%v\n", err)
			os.Exit(1)
		}
	}
}

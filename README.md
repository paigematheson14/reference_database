# reference_database
reference database manuscript code and notes



did this in the jupyter terminal environment 

```
ssh mahuika
```


# step 1: install biopython
this is tricky because of the nesi permissions. I used a roundabout way to install it using the following code which seemed to work:

```
python3 -m venv myenv
```
```
source myenv/bin/activate
```
```
pip install --upgrade pip
```
```
pip install biopython

```

# step 2: create a .py script to parse the .gb file into a csv file so that is more readable and easier to analyse

```
nano parse_genbank.py
```
```
from Bio import SeqIO
import csv

# Input and output files
genbank_file = "example.gb"
csv_file = "output.csv"

# Open CSV file for writing
with open(csv_file, mode="w", newline="") as csv_out:
    writer = csv.writer(csv_out)

    # Write CSV header
    writer.writerow([
        "Locus", "Sequence Length (bp)", "Molecule Type", "Topology", "Division", "Date",
        "Definition", "Accession", "Version", "BioProject", "BioSample", "Keywords",
        "Organism", "Taxonomy", "Reference Authors", "Reference Title", "Assembly Information",
        "Chromosome", "Sex", "Tissue Type", "Geographic Location", "Collection Date", "Breed"
    ])

    # Parse GenBank file
    for record in SeqIO.parse(genbank_file, "genbank"):
        # Extract core metadata
        locus = record.name
        seq_length = len(record.seq)
        molecule_type = record.annotations.get("molecule_type", "N/A")
        topology = record.annotations.get("topology", "N/A")
        division = record.annotations.get("data_file_division", "N/A")
        date = record.annotations.get("date", "N/A")
        definition = record.description
        accession = record.id
        version = record.annotations.get("sequence_version", "N/A")

        # Extract DBLinks (BioProject, BioSample)
        bioproject = biosample = "N/A"
        for dbxref in record.dbxrefs:
            if "BioProject" in dbxref:
                bioproject = dbxref.split(":")[1].strip()
            elif "BioSample" in dbxref:
                biosample = dbxref.split(":")[1].strip()

        # Extract keywords
        keywords = "; ".join(record.annotations.get("keywords", ["N/A"]))

        # Extract organism & taxonomy
        organism = record.annotations.get("organism", "N/A")
        taxonomy = "; ".join(record.annotations.get("taxonomy", []))

        # Extract reference authors and titles
        reference_authors = reference_title = "N/A"
        if record.annotations.get("references"):
            ref = record.annotations["references"][0]
            reference_authors = "; ".join(ref.authors.split(", ")) if ref.authors else "N/A"
            reference_title = ref.title if ref.title else "N/A"

        # Extract Assembly Information from the COMMENT section
        assembly_info = "N/A"
        comment_data = record.annotations.get("comment", "")
        if comment_data:
            comment_lines = comment_data.split("\n")
            assembly_lines = [line.strip() for line in comment_lines if "Assembly" in line]
            assembly_info = "; ".join(assembly_lines) if assembly_lines else "N/A"

        # Extract feature annotations
        chromosome = sex = tissue_type = geo_location = collection_date = breed = "N/A"
        for feature in record.features:
            if feature.type == "source":
                chromosome = feature.qualifiers.get("chromosome", ["N/A"])[0]
                sex = feature.qualifiers.get("sex", ["N/A"])[0]
                tissue_type = feature.qualifiers.get("tissue_type", ["N/A"])[0]
                geo_location = feature.qualifiers.get("geo_loc_name", ["N/A"])[0]
                collection_date = feature.qualifiers.get("collection_date", ["N/A"])[0]
                breed = feature.qualifiers.get("note", ["N/A"])[0] if "note" in feature.qualifiers else "N/A"

        # Write extracted data to CSV
        writer.writerow([
            locus, seq_length, molecule_type, topology, division, date,
            definition, accession, version, bioproject, biosample, keywords,
            organism, taxonomy, reference_authors, reference_title, assembly_info,
            chromosome, sex, tissue_type, geo_location, collection_date, breed
        ])

print(f"CSV file '{csv_file}' has been created successfully.")
```
```
python3 ./parse_genbank.py
```

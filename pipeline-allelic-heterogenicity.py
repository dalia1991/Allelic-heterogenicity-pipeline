from Bio import Entrez
import sys
import pickle
import MySQLdb
from ucsc.api import Genome, Sequence
import requests


def main():
    Entrez.email = "dalia.kamal.peds@gmail.com"
    disease = input("enter disease: ")
    with open("gene_list.txt","rb")as f :
        diseases = pickle.load(f)
    #adding disease to my list if not before 
    if disease not in diseases.keys():
        gene_list = get_gene_from_disease(disease)
        diseases[disease] = gene_list
        line = pickle.dumps(diseases)
        with open ("gene_list.txt","wb") as file:
            file.write(line)
    #for sake of demonstration i will use only gene_id_0 but can iteratea over all genes 
    gene_id_0 = diseases[disease][0]
    query_file = input("query fasta file: ")
    with open(query_file,"r") as f:
        lines = f.readlines()
        name,seq0=lines[0].split(">")
        seq = "".join(lines[1:])
        query_sequence = seq0 + seq
    for transcript in transcript_generator_for_alignment(gene_id_0):
        for trascript_variant in get_variants_in_transcript(transcript, genome="hg38"):
            chrom = transcript["chrom"]
            start = transcript["tx_start"]
            end = transcript["tx_end"]

            transcript_sequence = fetch_genome_sequence(chrom, start, end)
            transcript_sequence = ''.join(transcript_sequence.splitlines()[1:]).strip()

            print(f"\n🧬 BLASTing against transcript {transcript['transcript_id']}...")
            run_online_blast(query_sequence, transcript_sequence)
# retrieving genes related to disease
def get_gene_from_disease(disease):
    term = f"{disease} AND Homo sapiens[Organism]"
    handle = Entrez.esearch(db="gene", term=term)
    record = Entrez.read(handle)
    handle.close()

    gene_ids = record["IdList"]
    if not gene_ids:
        return None

    # Use esummary for quick info
    id_str = ",".join(gene_ids)
    handle = Entrez.esummary(db="gene", id=id_str)
    summary = Entrez.read(handle)
    handle.close()

    return [
        doc["Name"]
        for doc in summary["DocumentSummarySet"]["DocumentSummary"]
    ]


#retriving transcript coordinates of gene_id_0 as demonstration
#gene_name = diseases[disease][0]
def transcript_generator_for_alignment(gene_name, genome="hg38"):
    """this function takeas agene as input and output all possible transcript including normal isoforms due to alterntive splicing"""
    db = MySQLdb.connect(host="genome-mysql.soe.ucsc.edu",
                         user="genome", passwd="", db=genome)
    cursor = db.cursor()

    query = f"""
        SELECT name, chrom, strand, txStart, txEnd, exonStarts, exonEnds
        FROM refGene
        WHERE name2 = '{gene_name}'
        ORDER BY txStart;
    """
    cursor.execute(query)
    for tx_id, chrom, strand, tx_start, tx_end, exon_starts, exon_ends in cursor.fetchall():
        transcript_info = {
            "transcript_id": tx_id,
            "chrom": chrom,
            "strand": strand,
            "tx_start": tx_start,
            "tx_end": tx_end
        }

        yield transcript_info  # pause + return, resume next loop on next call
    db.close()

#selecting variants from column snps with clinical significance
def get_variants_in_transcript(transcript, genome="hg38"):
    """ this function retrieve the pathological variant of the isoform it takes an isoform as input and output its pathological variants"""
    db = MySQLdb.connect(host="genome-mysql.soe.ucsc.edu",
                         user="genome", passwd="", db=genome)
    cursor = db.cursor()

    query = f"""
        SELECT chrom, chromStart, chromEnd, name, refUCSC, func
        FROM snp151
        WHERE chrom = '{transcript["chrom"]}'
        AND chromStart >= {transcript["tx_start"]}
        AND chromEnd <= {transcript["tx_end"]}
        AND bitfields LIKE '%clinically-assoc%';
    """
    cursor.execute(query)
    variants = cursor.fetchall()
    
    variant_count = 0
    splice_site_count = 0
    transcript_variants = []

    for chrom, chromStart, chromEnd, name, refNCBI, func in variants:
        variant = {
            "chrom": chrom,
            "chromStart": chromStart,
            "chromEnd": chromEnd,
            "name": name,
            "refNCBI": refNCBI,
            "func": func
        }

        if "splice" in func.lower():
            splice_site_count += 1

        transcript_variants.append(variant)
        variant_count += 1

    db.close()
    
    print(f"ClinVar Variants: {variant_count} (including {splice_site_count} splice site mutation{'s' if splice_site_count != 1 else ''})")
    for variant in transcript_variants :
        yield variant

        

    
def fetch_genome_sequence(chrom, start, end, genome="hg38"):
    url = f"https://genome.ucsc.edu/cgi-bin/das/{genome}/dna?segment={chrom}:{start},{end}"
    response = requests.get(url)
    
    if response.status_code != 200:
        return f"Failed to fetch sequence: {response.status_code}"
    
    # Extract sequence from XML response
    lines = response.text.splitlines()
    seq_lines = [line.strip() for line in lines if not line.startswith("<")]
    sequence = "".join(seq_lines).upper()
    
    header = f">{genome}_{chrom}:{start}-{end}"
    return f"{header}\n{sequence}"



# using blast to align query and subject gene
import tempfile
from Bio.Blast import NCBIWWW, NCBIXML

def run_online_blast(query_sequence, subject_sequence, program="blastn", database="nt"):
    # Combine the subject sequence into the query payload
    # NCBI BLAST only allows querying a sequence against a database, not a specific sequence
    # So we treat the subject as the database (conceptually), but in practice — only query is used

    result_handle = NCBIWWW.qblast(program=program, database=database, sequence=query_sequence)

    blast_record = NCBIXML.read(result_handle)

    # Print summary of top alignment
    if blast_record.alignments:
        alignment = blast_record.alignments[0]
        hsp = alignment.hsps[0]
        print(f"\n🔍 Query vs Transcript")
        print(f"Transcript Approximate Match: {alignment.title}")
        print(f"E-value: {hsp.expect}")
        print(f"Identity: {hsp.identities}/{hsp.align_length}")
        print(f"\nQuery:\n{hsp.query}\nMatch:\n{hsp.match}\nSubject:\n{hsp.sbjct}")
    else:
        print("No hits found for transcript.")

    return blast_record

main()


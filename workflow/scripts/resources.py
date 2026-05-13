import os


class Resources:
    """Gets URLs and file names of fasta and GTF files for a given genome and build"""

    # create genome directory
    os.makedirs("resources/", exist_ok=True)

    def __init__(self, genome, build):
        self.genome = genome
        self.build = build

        # base URLs
        base_url = f"https://ftp.ensembl.org/pub/release-{build}/"

        if "hg" in genome:
            if genome == "hg19":
                name = "GRCh37"

                self.blacklist_url = "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/hg19-blacklist.v2.bed.gz"
            elif genome == "hg38":
                name = "GRCh38"

                self.blacklist_url = "https://www.encodeproject.org/files/ENCFF356LFX/@@download/ENCFF356LFX.bed.gz"

            # create URLs for genome files
            self.fasta_url = f"{base_url}fasta/homo_sapiens/dna/Homo_sapiens.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = (
                f"{base_url}gtf/homo_sapiens/Homo_sapiens.{name}.{build}.gtf.gz"
            )

        elif "mm" in genome:
            if genome == "mm38":
                name = "GRCm38"

                self.blacklist_url = "https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz"

            elif genome == "mm39":
                name = "GRCm39"
                self.blacklist_url = "https://github.com/niekwit/chip-seq/raw/refs/heads/main/resources/GRCm39_excluded.bed.gz"

            # create URLs for genome files
            self.fasta_url = f"{base_url}fasta/mus_musculus/dna/Mus_musculus.{name}.dna.primary_assembly.fa.gz"
            self.gtf_url = (
                f"{base_url}gtf/mus_musculus/Mus_musculus.{name}.{build}.gtf.gz"
            )

        elif "test" in genome:

            self.fasta_url = "https://github.com/niekwit/damid-seq/raw/main/.test_pe/Homo_sapiens.GRCh38.dna.primary_assembly_chr11.fa.gz"
            self.gtf_url = "https://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.gtf.gz"

        # downloaded unzipped file names
        self.fasta = self._file_from_url(self.fasta_url)
        self.gtf = self._file_from_url(self.gtf_url)
        self.blacklist = self._file_from_url(self.blacklist_url)

    def _file_from_url(self, url):
        """Returns file path for unzipped downloaded file"""

        return f"resources/{os.path.basename(url).replace('.gz','')}"

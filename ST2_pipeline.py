import os
import sys
import parseStringTieGtf

class ST2_pipeline:
    def __init__(self, args): #fq1, fq2, gtf_TE, bed_TE, genome_index, tissue, outFolder
        self.fq1, self.fq2, self.TEgtf, self.TEbed, self.genome_index, self.tissue, self.outFolder = args[1:]
        self.sam_big = self.outFolder + "/" + self.tissue + ".sam"
        self.bam_big = self.outFolder + "/" + self.tissue + ".bam"
        self.sorted_bam_big = self.outFolder + "/" + "sorted_" + self.tissue + ".bam"
        self.gtfST2 = self.outFolder + "/" + self.tissue + ".strtie2.gtf"
        self.abundancyST2 = self.outFolder + "/" + self.tissue + ".strtie2.abundancy"
        self.parsed_gtf = self.gtfST2 + "_modified_abundancyTab"
        self.hist_file = self.outFolder + "/" + self.tissue + ".hist"
        self.TEintsect_bam = self.outFolder + "/" + "TEintsect_" + self.tissue + ".bam"
        self.sortedTEintsect_bam = self.outFolder + "/" + "sorted_TEintsect_" + self.tissue + ".bam"
        self.unspliced_sorted_TEintsectHA = self.outFolder + "/" + "unspliced_sorted_TEintsect_" + self.tissue + ".bam"
        self.main()

    def main(self):
        ## 1. mapping
        print("1. Mapping")
        os.system("nohup /home/ikirov/Tools/hisat2-2.0.0-beta/hisat2 -k 200 -p 15 -x {0} -1 {1} -2 {2} -S {3}".format(
            self.genome_index, self.fq1, self.fq2, self.sam_big
        ))

        ## 2. sam to bam and sort
        print("2. Sam to Bam")
        os.system("samtools view -Sb {0} > {1}".format(self.sam_big, self.bam_big))
        os.system("rm {}".format(self.sam_big))
        print("Bam sorting")
        os.system("samtools sort -o {0} {1}".format(self.sorted_bam_big, self.bam_big))

        ## 3. StringTi2
        print("3. StringTie2 is running")
        os.system("/home/ikirov/Tools/stringtie2/stringtie {0} -u -o {1} -t -g 500 -p 15 -A {2} -v -G {3}".format(
            self.sorted_bam_big,self.gtfST2, self.abundancyST2, self.TEgtf
        ))

        ## 4. Get formatted gtf
        print("4. Parsing gtf file and {}_modified_abundancyTab generation".format(self.gtfST2))
        parseStringTieGtf.StringTIe2_parser(self.gtfST2)


        ## PREPARE for IGV
        print("5. INtersection for IGV visualisation and truncated bam generation")
        os.system("intersectBed -abam {0} -b {1} -wa -u -ubam > {2}".format(
            self.sorted_bam_big, self.TEbed, self.TEintsect_bam
        ))

        print("\t\t5.1Sorting intersect bam file")
        os.system("bamtools sort -in {0} -out {1}".format(self.TEintsect_bam, self.sortedTEintsect_bam))

        print("\t\t5.2Indexing intersect bam file")
        os.system("bamtools index -in {}".format(self.sortedTEintsect_bam))



        ##count coverage
        print("6. Coverage counting and {} file generation".format(self.hist_file))

        ## remove splice junctions
        print("\t\t6.1remove splice junctions...")
        os.system("nohup samtools view -h -F 4 {0} | awk '$6 !~ /N/ || $1 ~ /@/' | samtools view -bS - > {1}".format(
            self.sortedTEintsect_bam, self.unspliced_sorted_TEintsectHA
        ))

        print("\t\t6.2hist file generation...")
        os.system("coverageBed -b {0} -a {1} -hist > {2}".format(
            self.unspliced_sorted_TEintsectHA, self.TEbed, self.hist_file
        ))

        print("FINISH")



if __name__ == '__main__':
    ST2_pipeline(sys.argv)

"""
It takes gtf file and abundancy file from StringTie2
This script will create file with transcript IDs linked to the original TE loci

"""
from collections import defaultdict
class StringTIe2_parser:
    def __init__(self, gtf):
        self.gtf = gtf
        self.id_dic = {}
        self.header = "TE_id\tReference\tStart\tEnd\tGene_ID\tTranscript_id\tCoverage\tFPKM\tTPM\tisTE_exact\tNumber_exons\tTranscript_length\n"
        self.main()

    def __get_annotation_value(self, inp):
        # transcript_id "STRG.22.1"
        # print(inp.split(' "'))
        return inp.split(' "')[1][:-1]

    def parseGTF(self):
        """
        create dictionary transcript reference id (=TE id)) -> StringtIE2 ID al the rest from file
        :return:
        """

        with open(self.gtf) as infile, open(self.gtf + "_modified_abundancyTab", "w") as outFile:
            outFile.write(self.header)
            toWrtite = {}
            tr_length = defaultdict()
            tran_id_gene_id = {}
            exons = defaultdict() #transcript_id:number of exons
            is_TE_transcript = False
            cnt = 0
            for lines in infile:
                if not lines.startswith("#"):
                    sp = lines.rstrip().split("\t")
                    chromosome, start, end = sp[0], sp[3], sp[4]
                    annotation = sp[-1][:-1].rstrip().split(";")
                    gene_id = self.__get_annotation_value(annotation[0])
                    transcript_id = self.__get_annotation_value(annotation[1])

                    if sp[2] == "transcript":
                        cov, fpkm, tpm = [self.__get_annotation_value(i) for i in [annotation[-3], annotation[-2], annotation[-1]]]
                        if "reference_id" in sp[-1]:
                            #gene_id "STRG.22"; transcript_id "STRG.22.1"; reference_id "01|chr|HA412.bronze.20140814_989772_996600"; ref_gene_id "01|chr|HA412.bronze.20140814_989772_996600"; cov "0.002636"; FPKM "0.000581"; TPM "0.001236";
                            TE_id = self.__get_annotation_value(annotation[2])
                            self.id_dic[gene_id] = TE_id
                            is_TE_transcript = True
                        else:
                            is_TE_transcript = False
                        tran_id_gene_id[transcript_id] = [gene_id, "\t".join([chromosome, start, end, gene_id, transcript_id, cov, fpkm, tpm, str(is_TE_transcript)])]

                    if sp[2] == "exon":
                        if transcript_id in exons:
                            exons[transcript_id] += 1
                            tr_length[transcript_id] += abs(int(start) - int(end))
                        else:
                            exons[transcript_id] = 1
                            tr_length[transcript_id] = abs(int(start) - int(end))

            for transcripts in tran_id_gene_id:
                if tran_id_gene_id[transcripts][0] in self.id_dic:
                    outFile.write(self.id_dic[tran_id_gene_id[transcripts][0]] + "\t"  +
                                  tran_id_gene_id[transcripts][1] + "\t" + str(exons[transcripts]) + "\t" + str(tr_length[transcripts]) + "\n")

                    cnt += 1


            print(cnt, "transcripts of TE have been collected")

            ### Add length column !!!!!!!!!!!!!!!

    def main(self):
        self.parseGTF()


#StringTIe2_parser(r"D:\Programming\R_working\retrotranscriptome2019\StringTie2_files\HA_matleaf_strtie2.gtf")
#StringTIe2_parser(r"D:\Programming\R_working\retrotranscriptome2019\StringTie2_files\HA_ligule_strtie2.gtf")
#StringTIe2_parser(r"D:\Programming\R_working\retrotranscriptome2019\StringTie2_files\HA_roots_strtie2.gtf")
#
# import os
#
# root_folder = r'D:\\VNIISB_2019\\retrotranscriptome\\Sunflower\\StringTie2_files\\'
# for files in os.listdir(root_folder):
#     if files.endswith("hist"):
#         print(root_folder + files)
#         #StringTIe2_parser(root_folder + files)
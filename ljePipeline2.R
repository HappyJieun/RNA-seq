#############################
# jieun Lee
# RNA-seq pipeline
# 2021.12.27 ~
#############################

setRepositories(ind = (1:8))

library(stringr)
library(progress)
library(tictoc)

# get project folder name
project <- commandArgs(trailingOnly = T)
#project <- "PRJNA494727"

print("***This pipeline is RNA sequencing for mouse(mus_musculus).***")
print(paste0("- Project : ", project))

# directory
dir_RNAseq <- "/disk3/bilje/RNA-seq"

setwd(dir_RNAseq)

dir_sample <- paste0(dir_RNAseq, "/0.Sample/", project)
dir_trim <- paste0(dir_RNAseq, "/1.Trimmometic/", project)
dir_ref <- paste0(dir_RNAseq, "/2.Reference/", project)
dir_SAM <- paste0(dir_RNAseq, "/3.SAM/", project)
dir_BAM <- paste0(dir_RNAseq, "/4.BAM/", project)
dir_sortBAM <- paste0(dir_RNAseq, "/5.Sorted_BAM/", project)
dir_quan <- paste0(dir_RNAseq, "/6.Quantification/", project)

# make project directory
system(paste0("mkdir ", dir_trim))
system(paste0("mkdir ", dir_trim, "/log"))
system(paste0("mkdir ", dir_ref))
system(paste0("mkdir ", dir_ref, "/log"))
system(paste0("mkdir ", dir_SAM))
system(paste0("mkdir ", dir_SAM, "/log"))
system(paste0("mkdir ", dir_BAM))
system(paste0("mkdir ", dir_BAM, "/log"))
system(paste0("mkdir ", dir_sortBAM))
system(paste0("mkdir ", dir_sortBAM, "/log"))
system(paste0("mkdir ", dir_quan))
system(paste0("mkdir ", dir_quan, "/log"))

#save sample's names
fastq <- list.files(dir_sample, pattern = ".fastq.gz")
fastq_name <- c()
for (i in 1:length(fastq)){
  fastq_name[i] <- strsplit(fastq[i], "_")[[1]][1]
}

# check number of sample 
fastq_list <- unique(fastq_name)
if (!((length(fastq_name))/2==length(fastq_list))){
  print("The samples are not paired-end reades.")
  quit()
}

#trimmometic
setwd(dir_trim)
print(paste0("1. Trimmometic of ", project, " project."))
tic("[1] \"trimmometic\"")
pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                     max = length(fastq_list), # Maximum value of the progress bar
                     style = 3,    # Progress bar style (also available style = 1 and style = 2)
                     width = 50,   # Progress bar width. Defaults to getOption("width")
                     char = "=")
trim <- c()
for (i in 1:length(fastq_list)){
  trim[i] <- paste0("java -jar /program/Trimmomatic/trimmomatic-0.39.jar PE -threads 16 -phred33 ", 
                    dir_sample, "/", fastq_list[i], "_1.fastq.gz ", dir_sample, "/", fastq_list[i], "_2.fastq.gz ",
                    dir_trim, "/", fastq_list[i], "_1P.fastq.gz ", dir_trim, "/", fastq_list[i], "_1uP.fastq.gz ", 
                    dir_trim, "/", fastq_list[i], "_2P.fastq.gz ", dir_trim, "/", fastq_list[i], "_2uP.fastq.gz ",
                    "ILLUMINACLIP:/program/Trimmomatic/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 2> ", dir_trim, "/log/", fastq_list[i], "_trim.log")
  system(trim[i])
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
  # system(paste0("Finished ", fastq_list[i], " trimming."))
}
close(pb) # Close the connection
toc()

print("Complete")

#reference & annotation
setwd(dir_ref)

speices1 <- "mus_musculus"
speices2 <- str_to_title(speices1)
ID <- "GRCm39"

#reference genome
print("2. Downloading reference genome file of mouse")
ref_gene <- paste0("wget 2> ", dir_ref, "/log/GRCm39_reference.log ",
                   "http://ftp.ensembl.org/pub/release-105/fasta/", speices1,
                   "/dna/", speices2, ".", ID, ".dna.primary_assembly.fa.gz")
tic("[1] \"Down reference genome\"")
system(ref_gene)
toc()
print("Reference genome file is saved")

#genome annotation
print("3. Downloading genome annotaion file of mouse.")
annote_gene <- paste0("wget 2> ", dir_ref, "/log/GRCm39_annotaion.log ",
                      "http://ftp.ensembl.org/pub/release-105/gtf/", speices1, "/", speices2, ".", ID, ".105.gtf.gz")
tic("[1] \"Down genome annotaion\"")
system(annote_gene)
toc()
print("Genome annotaion file is saved")

# unzip
setwd(dir_ref)
print("4. Unzip")
tic("[1] \"unzip\"")
system("gzip -d *gz")
toc()
print("Complete")

#reference genome indexing
print("5. Indexing reference genome")
ref_index <- paste0("/program/HISAT2/hisat2-build 2> ", dir_ref, "/log/", ID, "_index.log ", 
                    dir_ref, "/", speices2, ".", ID, ".dna.primary_assembly.fa ", dir_ref, "/", speices2, "_ind")
tic("[1] \"Indexing\"")
system(ref_index)
toc()
print("Complete")

#fastq to SAM
setwd(dir_SAM)
print("6. fastq to SAM")
pb <- txtProgressBar(min = 0,     
                     max = length(fastq_list),
                     style = 3,   
                     width = 50,
                     char = "=")
SAM <- c()
tic("[1] \"SAM\"")
for (i in 1:length(fastq_list)){
  SAM[i] <- paste0("/program/HISAT2/hisat2 2> ", dir_SAM, "/log/", fastq_list[i], "_sam.log ", 
                   "-x ", dir_ref, "/", speices2, "_ind ",
                   "-1 ", dir_trim, "/", fastq_list[i], "_1P.fastq.gz ",
                   "-2 ", dir_trim, "/", fastq_list[i], "_2P.fastq.gz ",
                   "-S ", dir_SAM, "/", fastq_list[i], ".sam")
  system(SAM[i])
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)
toc()
print("Complete")

#SAM to BAM
setwd(dir_BAM)
print("7. SAM to BAM")
pb <- txtProgressBar(min = 0,     
                     max = length(fastq_list),
                     style = 3,   
                     width = 50,
                     char = "=")
BAM <- c()
tic("[1] \"BAM\"")
for (i in 1:length(fastq_list)){
  BAM[i] <- paste0("/program/samtools/bin/samtools view -bS ",
                   dir_SAM, "/", fastq_list[i], ".sam",
                   " > ", dir_BAM, "/", fastq_list[i], ".bam 2> ", dir_BAM, "/log/", fastq_list[i], "_bam.log")

  system(BAM[i])
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)
toc()
print("Complete")

#BAM to Soreted BAM
setwd(dir_sortBAM)
print("8. BAM to Soreted BAM")
pb <- txtProgressBar(min = 0,     
                     max = length(fastq_list),
                     style = 3,   
                     width = 50,
                     char = "=")
sorted_BAM <- c()
tic("[1] \"Soreted BAM\"")
for (i in 1:length(fastq_list)){
  sorted_BAM[i] <- paste0("/program/samtools/bin/samtools sort ", dir_BAM, "/", fastq_list[i], ".bam ",
                          "-o ", dir_sortBAM, "/", fastq_list[i], "_sorted.bam 2> ", dir_sortBAM, "/log/", fastq_list[i], "_sorted.log")
  system(sorted_BAM[i])
  Sys.sleep(0.1)
  setTxtProgressBar(pb, i)
}
close(pb)
toc()
print("Complete")

#Quantification
setwd(dir_quan)
print("9. Quantification")
pb <- txtProgressBar(min = 0,     
                     max = length(fastq_list),
                     style = 3,   
                     width = 50,
                     char = "=")
feature_count <- c()
tic("[1] \"FeatureCount\"")
for (i in 1:length(fastq_list)){
   feature_count[i] <- paste0("/program/subread/bin/featureCounts -p -a ", dir_ref, "/", speices2, ".", ID, ".105.gtf ",
                             "-t exon -g gene_id -s 0 ", "-o ",
                             dir_quan, "/", fastq_list[i], "_fc.txt ", dir_sortBAM, "/", fastq_list[i], "_sorted.bam ",
                             "2> ", dir_quan, "/log/", fastq_list[i], "_fc.log")
   system(feature_count[i])
   Sys.sleep(0.1)
   setTxtProgressBar(pb, i)
}
close(pb)
toc()

print("Complete")
print("Well done!")

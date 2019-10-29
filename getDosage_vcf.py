import gzip

for i in range(1, 23):
    outName = "/project2/yangili1/zpmu/dice/data/fastqtl/chr%i.dose.r2.nodup.rename.filtered.maf5e-2.recode.dosage.gz"%(i)
    inName = "/project2/yangili1/zpmu/dice/data/fastqtl/chr%i.dose.r2.nodup.rename.filtered.maf5e-2.recode.recode.vcf.gz"%(i)
    
    with gzip.open(inName, "rt") as inFile, gzip.open(outName, "wt") as outFile:
        for ln in inFile:
            if ln.startswith("##"):
                continue # skip annotation lines
            if ln.startswith("#CHROM"):
                myVec = ln.rstrip().split("\t")
                iid = myVec[9:]
                newLn = "\t".join(iid)
                outFile.write("sid\tALT\t" + newLn + "\n")
                continue
            myVec = ln.rstrip().split("\t")
            sid = myVec[2]
            ref = myVec[4]
            count = [x.split(":")[1] for x in myVec[9:]]
            newVec = [sid] + [ref] + count
            newLn = "\t".join(newVec)
            outFile.write(newLn + "\n")


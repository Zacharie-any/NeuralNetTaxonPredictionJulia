using BioSequences, Flux, DataFrames, BSON


include("MyBioTool.jl")

primer341F = rna"CCUACGGGNGGCWGCAG"
primer785R = rna"GACUACHVGGGUAUCUAAUCC" #example primers (here targeting v3-v4 region of 16S rrna)

taxonnames=["Firmicutes","Proteobacteria","Cyanobacteria","Metazoa"] #example

T,taxid2name = import_ref("data/SILVA_138.1_SSURef_NR99_tax_silva.fasta","data/SILVA_138.1_SSURef_Nr99.full_metadata","data/tax_slv_ssu_138.1.tre","data/tax_slv_ssu_138.1.acc_taxid","data/tax_slv_ssu_138.1.map")
println("data imported successfuly")

T = T[1001:2000,:]

labeling!(T;taxonnamespool=taxonnames,taxid2name=taxid2name)
println("data labeled")
dropmissing!(T,:label)

T = harmonize_size_by_label(T;minnumelem=2)
println("data harmonized")  
println("$(nrow(T)) records") 

p = pre_process_learning_sets(T;splitpercentage=0.9,wordsize=3,samplesize=300, coverage=5)
println("data preprocessed successfuly")
#T=nothing #freeing memory

BSON.@save "preprocessed/traintest.BSON" p
(Xtrain,Ytrain,Xtest,Ytest,wordencoding,trainoffsets,testoffsets) = p


model = learn(Xtrain,Ytrain,Xtest,Ytest;numconvlayer=2,filtersize=9,filternum=32,embeddingdims=64,maxepoch=2,batchsize=128,savepath = "runs/1/",savename="model.bson",encoding=wordencoding)
println("model completed the learning step")


using BioSequences, Flux, DataFrames, Plots


include("MyBioTool.jl")

# plot the error rate with the offset : 
#   set 0.3 splitpercentage
#   set a big coverage
#   set in shotgun mode

# shotgun mode vs with primers : memory usage also

# plots accuracy with number of params (filtersize ->| and numconvlayer ->|)

# plots accuracy with samplesize

# plots accuracy with noisepercentage

# plots accuracy with taxon distance 
#   in node distance
#   in handpicked names for special learning

Ts,taxid2name = import_ref("data/SILVA_138.1_SSURef_NR99_tax_silva.fasta","data/SILVA_138.1_SSURef_Nr99.full_metadata","data/tax_slv_ssu_138.1.tre","data/tax_slv_ssu_138.1.acc_taxid","data/tax_slv_ssu_138.1.map")
println("data imported successfuly")

taxonnames=["Firmicutes","Proteobacteria","Cyanobacteria","Metazoa"] # ! Metazoa has no v3v4 primer region by the primer so it will be let go in preprocessing

    #1 

    T=Ts
    labeling!(T;taxonnamespool=taxonnames,taxid2name=taxid2name)
    println("data labeled")
    dropmissing!(T,:label)

    T = harmonize_size_by_label(T;minnumelem=2)
    println("data harmonized")
    println("$(nrow(T)) records")

    Xtrain,Ytrain,Xtest,Ytest,wordencoding,trainoffsets,testoffsets = pre_process_learning_sets(T;splitpercentage=0.3,samplesize=300,coverage=50,wordsize=3)
    println("data preprocessed successfuly")

    model = learn(Xtrain,Ytrain,Xtest,Ytest;filtersize=9,embeddingdims=64,maxepoch=10,batchsize=32,encoding=wordencoding,does_shuffle = false)
    println("model completed the learning step")

    primer341F = rna"CCUACGGGNGGCWGCAG"
    primer785R = rna"GACUACHVGGGUAUCUAAUCC"

    plot_accuracy_v_offset(model,Xtest,Ytest,testoffsets,20)
    re,r1,r2=get_primersregion(T.Sequence,primer341F,primer785R)
    plot!(r1:r2,[1 for i in r1:r2],label = "v3-v4 region")

#2

    


    labeling!(T;nodedistancetodomain=2)
    println("data labeled")

    harmonize_size_by_label!(T;minnumelem=2)
    println("data harmonized")

    Xtrain,Ytrain,Xtest,Ytest,wordencoding,trainoffsets,testoffsets = pre_process_learning_sets(T;splitpercentage=0.5,wordsize=3,onprimers=[primer341F,primer785R])
    println("data preprocessed successfuly")

    model = learn(Xtrain,Ytrain;numconvlayer=3,filtersize=12,embeddingdims=64,maxepoch=40,batchsize=false,encoding=wordencoding)
    println("model completed the learning step")
    acc2 = measure_accuracy(model,Xtest,Ytest)

#3

    labeling!(T;nodedistancetodomain=2)
    println("data labeled")

    harmonize_size_by_label!(T;minnumelem=2)
    println("data harmonized")

    Xtrain,Ytrain,Xtest,Ytest,wordencoding,trainoffsets,testoffsets = pre_process_learning_sets(T;splitpercentage=0.5,samplesize=300,coverage=3,wordsize=3)
    println("data preprocessed successfuly")
    acc3=Matrix([])
    for filtersize in [3,6,9,12] 
    for numconvlayer in [2,3,4]
    model = learn(Xtrain,Ytrain;numconvlayer,filtersize,embeddingdims=64,maxepoch=40,batchsize=false,encoding=wordencoding)
    println("model completed the learning step")

    end
    end

#4


    T=Ts
    labeling!(T;taxonnamespool=taxonnames,taxid2name=taxid2name)
    println("data labeled")
    dropmissing!(T,:label)

    T = harmonize_size_by_label(T;minnumelem=2)
    println("data harmonized")
    println("$(nrow(T)) records")

    
    acc4=[]
    optacc4=[]
    for samplesize in [50,100,200,300] 
    
        Xtrain,Ytrain,Xtest,Ytest,wordencoding,trainoffsets,testoffsets = pre_process_learning_sets(T;splitpercentage=0.3,samplesize=samplesize,coverage=50,wordsize=3)
        println("data preprocessed successfuly")
    
        model = learn(Xtrain,Ytrain,Xtest,Ytest;filtersize=9,embeddingdims=64,maxepoch=10,batchsize=32,encoding=wordencoding,does_shuffle = false,savename="savename = model$(samplesize).bson")
        println("model completed the learning step")
        
        append!(acc4,measure_accuracy(model,Xtest,Ytest))
        #append!(optacc4,measure_bayes_optimal_accuracy(Xtrain,Ytrain))
    end
    plot([50,100,200,300] ,acc4,xlabel="samplesize",ylabel="accuracy")

#5
    acc5=[]
    for nodedistancetodomain in [0,1,2,3,4] 
    labeling!(T;nodedistancetodomain)
    println("data labeled")

    harmonize_size_by_label!(T;minnumelem=2)
    println("data harmonized")
    

    Xtrain,Ytrain,Xtest,Ytest,wordencoding,trainoffsets,testoffsets = pre_process_learning_sets(T;splitpercentage=0.5,samplesize=300,coverage=3,wordsize=3)
    println("data preprocessed successfuly")
    

    model = learn(Xtrain,Ytrain;numconvlayer,filtersize,embeddingdims=64,maxepoch=40,batchsize=false,encoding=wordencoding)
    println("model completed the learning step")
    append!(acc5,measure_accuracy(model,Xtest,Ytest))
    
    end

#6

    acc6=[]
    optacc6=[]
    for noisepercentage in [0,0.001,0.01,0.05] 
    T=Ts
    labeling!(T;taxonnamespool=taxonnames,taxid2name=taxid2name)
    println("data labeled")
    dropmissing!(T,:label)

    T = harmonize_size_by_label(T;minnumelem=2)
    println("data harmonized")
    println("$(nrow(T)) records")

    introduce_noise!(T;noisepercentage)

    Xtrain,Ytrain,Xtest,Ytest,wordencoding,trainoffsets,testoffsets = pre_process_learning_sets(T;splitpercentage=0.9,samplesize=300,coverage=3,wordsize=3)
    println("data preprocessed successfuly")

    model = learn(Xtrain,Ytrain;numconvlayer,filtersize,embeddingdims=64,maxepoch=3,batchsize=false,encoding=wordencoding)
    println("model completed the learning step")
    append!(acc6,measure_accuracy(model,Xtest,Ytest))
    append!(optacc6,measure_bayes_optimal_accuracy(Xtrain,Ytrain))
    end

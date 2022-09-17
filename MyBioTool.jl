using FASTX
using BioSequences
using Flux
using DataFrames
using Random
using CSV
using Plots
import BSON
using ProgressMeter: @showprogress

function import_ref(fastapath,metapath,treepath,acctaxidpath,taxidnamepath)
    """
    imports relevent reference sequences and taxonomic information from SILVA SSU 138.1 files
    """
    
    function get_table(fastapath)
        s = []
        a = []
        open(FASTA.Reader, fastapath) do reader
            for record in reader
                append!(a,[FASTA.identifier(record)])
                append!(s,[FASTA.sequence(record)])
            end
        end
        return DataFrame(Accession_number = a, Sequence = s)
    end

    function filter_16S!(T,metapath)
       
        metadata = DataFrame(CSV.File(metapath))
        Acc = []
        for i in 1:nrow(metadata)
            append!(Acc,[(x->join([x.acc,x.start,x.stop],"."))(metadata[i,["acc","start","stop"]])])
        end

        metadata[!,"Accession_number"] = Acc
        return filter!(:product =>  x -> if ~ismissing(x) occursin("16S", x) else false end,innerjoin(T,metadata, on =:Accession_number))[!,["Accession_number","Sequence"]]
    end

    function get_parents(taxid,treepath)
        function get_father(taxid,char_tree)
            

            k = match(Regex(join(["[(),]",string(taxid),"[(),]"])),char_tree).offset
            k = k+1
            opened_parentheses = 0
            while ~(opened_parentheses == 0 && char_tree[k] == ')')
                if char_tree[k] == '('
                    opened_parentheses = opened_parentheses + 1
                end
                if char_tree[k] == ')'
                    opened_parentheses = opened_parentheses - 1
                end
                k = k + 1
            end
            k = k + 1
            start = k
            while ~(char_tree[k] == ')' || char_tree[k] == '(' || char_tree[k] == ','|| char_tree[k] == ';')
                k = k + 1;
            end
            stop = k-1

            return parse(Int,char_tree[start:stop])

        end
        parents = []
        open(treepath)do file 
            char_tree = read(file,String)
            

            node = taxid
            
            while ~(node in [2,3,4])
                append!(parents,node)
                node = get_father(node,char_tree)
            end
            append!(parents,node)
            
            
        end
        return parents
    end

    

    function get_taxid!(T,acctaxidpath)
        acc2taxid = Dict()
        open(acctaxidpath) do file
            for ln in eachline(file)
                a,b=split(ln,"\t")
                acc2taxid[a]=parse(Int,b)
            end
        end
        T[!,"taxid"] = map(T[!,"Accession_number"])do x return (acc2taxid[x])end
    end

    function get_taxnames!(T,taxidnamepath)
        taxid2name = Dict()
        open(taxidnamepath) do file
            for ln in eachline(file)
                a,b=split(ln,"\t")
                taxid2name[parse(Int,a)]=b
            end
        end
        T[!,"name"] = map(T[!,"taxid"])do x return (taxid2name[x])end
        return taxid2name
    end
    
    T = get_table(fastapath)
    println("fasta file read. $(nrow(T)) records.")

    

    T = filter_16S!(T,metapath)
    println("16S sequences only filtered. $(nrow(T)) records.")
    
    get_taxid!(T,acctaxidpath)
    println("taxid added")
    
    T[!,"parents"] = @showprogress map(T[!,"taxid"])do x 
        return (get_parents(x,treepath)) end
    println("parents added")
    
    taxid2name = get_taxnames!(T,taxidnamepath)
    println("scientific name added")

    

    return T, taxid2name

end

function make_refseq!(T::AbstractDataFrame)
    println(eltype(typeof(T[1,"Sequence"])))
    if eltype(typeof(T[1,"Sequence"])) == DNA
        base = [DNA_A,DNA_C,DNA_G,DNA_T,DNA_N] 
        null = DNA_N
    elseif eltype(typeof(T[1,"Sequence"])) == RNA
        base = [RNA_A,RNA_C,RNA_G,RNA_U,RNA_N] 
        null = RNA_N
    end
    @showprogress for i in 1:length(T[!,"Sequence"])
        for j in 1:length(T[i,"Sequence"])
            if ~(T[i,"Sequence"][j] in base)
                T[i,"Sequence"][j] = null
            end
        end
        # T[i,"Sequence"]=BioSequences.ReferenceSequence(T[i,"Sequence"])
    end
end
function make_refseq(Sequences)
    println(eltype(typeof(Sequences[1])))
    if eltype(typeof(Sequences[1])) == DNA
        base = [DNA_A,DNA_C,DNA_G,DNA_T,DNA_N] 
        null = DNA_N
    elseif eltype(typeof(Sequences[1])) == RNA
        base = [RNA_A,RNA_C,RNA_G,RNA_U,RNA_N] 
        null = RNA_N
    end
    @showprogress for i in eachindex(Sequences)
        for j in eachindex(Sequences)
            if ~(Sequences[i][j] in base)
                Sequences[i][j] = null
            end
        end
    end
    return Sequences
end

function get_primersregion(Sequences,primerA,primerB)
    "returns (if it exists) the sequences comprised between the primers (direct and reverse) "
    println("primer A :$primerA , $(typeof(primerA))")
    println("primer B :$primerB , $(typeof(primerA))")
    function find_primers_region(sequence)
        

        query1 = ApproximateSearchQuery(primerA, iscompatible)
        range1 = findfirst(query1, 0, sequence)

        if !(isnothing(range1))
            query2 = ApproximateSearchQuery(primerB, iscompatible)
            range2 = findfirst(query2, 0, BioSequences.reverse_complement(sequence))
            if !(isnothing(range2))
                return sequence[range1[1]:(length(sequence)-range2[1]+1)],range1[1],length(sequence)-range2[1]
            end
        end
        query1 = ApproximateSearchQuery(primerB, iscompatible)
        range1 = findfirst(query1, 0, sequence)

        if !(isnothing(range1))
            query2 = ApproximateSearchQuery(primerA, iscompatible)
            range2 = findfirst(query2, 0, BioSequences.reverse_complement(sequence))
            if !(isnothing(range2))
                return sequence[range1[1]:(length(sequence)-range2[1])],range1[1],length(sequence)-range2[1]
            end
        end

        return nothing, nothing, nothing

    end
    
    v3v4s = []
    local rangg1
    local rangg2
    for seq in Sequences
        res,r1,r2 = find_primers_region(seq)
        if !(isnothing(res))
            append!(v3v4s,[res])
            rangg1,rangg2=r1,r2

        else 
            append!(v3v4s,[missing])
        end
    end
    return v3v4s,rangg1,rangg2
end

function labeling!(T;taxonnamespool=false,taxonidspool=false,nodedistancetodomain=false,taxid2name=false)
    """
    adds a label field to T which will be used during the learning step.
    """
    if ~isa(taxonnamespool,Bool) && ~isa(taxid2name,Bool)

        T[!,"label"] = map(T[!,"parents"])do parents 
            for node in parents 
                if taxid2name[node] in taxonnamespool
                    return node
                else
                    continue
                end
            end
            return missing
        end
    elseif ~isa(taxonidspool,Bool)
        T[!,"label"] = map(T[!,"parents"])do parents 
            for node in parents 
                if node in taxonidspool
                    return node
                else
                    continue
                end
            end
            return missing
        end

    elseif ~isa(nodedistancetodomain,Bool)
        T[!,"label"] = map(T[!,"parents"])do parents 
            if length(parents)>=length(parents)-nodedistancetodomain>0
                return parents[length(parents)-nodedistancetodomain]
            else
                return missing
            end
        end

    end
end

function harmonize_size_by_label(T;minnumelem=2)
    df = combine(groupby(T,:label),nrow => :count)
    selectedlabels = df[df.count .>= minnumelem,:label]

    count_=Dict()
    for lab in unique(T[!,"label"])
        count_[lab]=0
    end

    return T[map(T[!,"label"])do lab begin 
        count_[lab]+=1
        
        return ((lab in selectedlabels) && count_[lab] <= minnumelem )
        
    end end,:]
end

function introduce_noise!(T;noisepercentage=0.001)
    if eltype(typeof(T[1,"Sequence"])) == DNA
        base = [DNA_A,DNA_C,DNA_G,DNA_T,DNA_N] 
    elseif eltype(typeof(T[1,"Sequence"])) == RNA
        base = [RNA_A,RNA_C,RNA_G,RNA_U,RNA_N] 
    end
    for i in 1:nrow(T)
        for j in 1:length(T[i,"Sequence"])
            if rand() < noisepercentage
                T[i,"Sequence"][j]=rand(base)
            end
        end
    end
end

function pre_process_learning_sets(T;splitpercentage=0.9,samplesize=300,coverage=3,wordsize=3,onprimers=false,initnumsamples=convert(Int,round(coverage*nrow(T)*length(T[1,"Sequence"])/samplesize*1.5)))
    
    coding_dict=Dict()
    N=1

    function splitting_sequences_on_primers()
        T[!,"SequencePrimersRegion"],r1,r2 = get_primersregion(T[!,"Sequence"],onprimers[1],onprimers[2])
        oldnrow = nrow(T)
        newnrow = nrow(dropmissing!(T,"SequencePrimersRegion"))
        println("Found $newnrow on $oldnrow samples containing a region within the given primers")
        primersregion = T[!,"SequencePrimersRegion"]
        primersregion = make_refseq(primersregion)
        println("converted to ReferenceSequence.")
        
        if eltype(primersregion[1]) == DNA
            null = dna"N"
        elseif eltype(primersregion[1]) == RNA
            null = rna"N"
        end
        samplesize = maximum(map(primersregion)do x length(x) end)
        X,Y = Matrix{Int}(undef,nrow(T),samplesize-wordsize+1),Vector{Int}(undef,nrow(T))
        offsets = Vector{Int}(undef,nrow(T))
        sn=1
        @showprogress for i in 1:nrow(T)
            sequence = primersregion[i]
            X[sn,:]=encode(sequence,null)
            Y[sn]=T[i,"label"]
            sn+=1
        end
        return X,Y,offsets
    end


    function splitting_sequences()

        make_refseq!(T)
        println("converted to ReferenceSequence.")

        X,Y = Matrix{Int}(undef,initnumsamples,samplesize-wordsize+1),Vector{Int}(undef,initnumsamples)
        offsets = Vector{Int}(undef,initnumsamples)
        sn=1

        @showprogress for i in 1:nrow(T)
            
            sequence = T[i,"Sequence"]
            S = 0
            while S < coverage * length(sequence)
                # if isa(samplesize,AbstractArray)
                #     sampsize = rand 
                # we can try to use random sample size and fill on the sides with NNN... seqs
                n = rand(1:(length(sequence)-samplesize+1))
                X[sn,:]=encode(sequence[n:n+samplesize-1],nothing)
                Y[sn]=T[i,"label"]
                offsets[sn]=n
                sn += 1
                S += samplesize
            end
            
        end

        X=X[1:sn-1,:]
        Y=Y[1:sn-1]
        offsets=offsets[1:sn-1]

        return X,Y,offsets
    end

    function encode(sequence,null)
        if length(sequence) < samplesize
            while length(sequence) != samplesize
                append!(sequence,null)
            end
        end
        coded_sequence = Vector{Int}(undef,samplesize-wordsize+1)
        
        for i in 1:(samplesize-wordsize+1)
            word = sequence[i:i+wordsize-1]
            if !(word in keys(coding_dict))
                coding_dict[word] = N
                N +=1
            end
            coded_sequence[i]=coding_dict[word]
        end
        return coded_sequence
    end

    function display_label_repartition(Y)
        

        df = DataFrame(label=Y)
        gdf = groupby(df,:label)
        labelcount = combine(gdf, nrow => :count)[!,:count]

        @info "label repartition in Ytrain :" sort(labelcount/length(Y))
    end

    function maketraintest(X,Y,offsets)
        function splitdf(df, pct)
            @assert 0 <= pct <= 1
            ids = collect(axes(df, 1))

            sel=ids.==0
            
            for type in unique(Y)
                
                typids = ids[df.y.==type]
                shuffle!(typids)
                tsel = collect(1:length(typids)).<=length(typids).*pct
                trtypids = typids[tsel]
                
                for t in trtypids
                    
                    sel = sel .| (ids.==t)
                end
                
            end
            
            
            return X[sel,:],X[.!sel,:], df[ sel, :],df[ .!sel, :]
        end
        df = DataFrame(y=Y,off=offsets)
        Xtrain,Xtest,dftr, dfte= splitdf(df, splitpercentage)
        Ytrain,trainoffsets,Ytest,testoffsets = dftr[!,"y"],dftr[!,"off"],dfte[!,"y"],dfte[!,"off"]
        X,Y,offsets,=nothing,nothing,nothing
        

        return  Xtrain,Ytrain,Xtest,Ytest,trainoffsets,testoffsets
    end

    

    T = T[randperm(nrow(T)),:]
    println("splitting and encoding step")

    if isa(onprimers,AbstractArray)
        @time X,Y,offsets = splitting_sequences_on_primers()
    else
        @time X,Y,offsets = splitting_sequences()
    end
    
    
    # T=nothing #free memory
    Xtrain,Ytrain,Xtest,Ytest,trainoffsets,testoffsets = maketraintest(X,Y,offsets)

    
    
    println("label repartition : ")
    display_label_repartition(Ytrain)

    return Xtrain,Ytrain,Xtest,Ytest,coding_dict,trainoffsets,testoffsets
end

function learn(Xtrain,Ytrain,Xtest,Ytest; numconvlayer=2, filtersize=9, filternum=32, poolsize=12, embeddingdims=64, maxepoch=20, batchsize=convert(Int,round(size(Xtrain)[1]/10)), encoding, η = 1e-3 ,λ = 1e-4,infotime = 2, savepath = "runs/",savename="model.bson", checktime = 5, does_shuffle=false)

    # η            learning rate  smaller -> more precise, slower
    # λ                  L2 regularizer param, implemented as weight decay
    # batchsize       batch size
    # epochs          number of epochs
    # infotime  	     report every `infotime` epochs
    # checktime      Save the model every `checktime` epochs. Set to 0 for no checkpoints. 
    # savepath    results path

    

    function report(epoch)
        train = eval_loss_accuracy(train_loader, model)
        test = eval_loss_accuracy(test_loader, model)        
        println("Epoch: $epoch   Train: $(train)   Test: $(test)")
        
    end
    
    @assert length(unique(Ytrain))==length(unique(Ytest))
    nclasses = length(unique(Ytrain))
    

    ytrain = Flux.onehotbatch(Ytrain, unique(Ytrain))
    ytest = Flux.onehotbatch(Ytest,unique(Ytest))

    Xtest = permutedims(Xtest)
    Xtrain = permutedims(Xtrain)

    samplesize = size(Xtrain)[1]

    Xtrain = reshape(Xtrain, size(Xtrain)[1], 1, size(Xtrain)[2])
    Xtest = reshape(Xtest, size(Xtest)[1], 1, size(Xtest)[2])

    train_loader = Flux.Data.DataLoader((Xtrain, ytrain), batchsize=batchsize, shuffle=does_shuffle)
    test_loader = Flux.Data.DataLoader((Xtest, ytest), batchsize=batchsize, shuffle=does_shuffle)

    layers = []
    append!(layers,[
        Flux.Embedding(length(encoding)=>embeddingdims),  
    ])

    append!(layers,[
            Conv((embeddingdims,filtersize),1=>32,relu,pad=SamePad())
            x->Flux.normalise(x,dims=3)
            Conv((1,filtersize),32=>64,relu,pad=SamePad())
            x->Flux.normalise(x,dims=3)
        ])
   

    # for k in 1:numconvlayer
    #     append!(layers,[
    #         Conv(((k==1 ? embeddingdims : 1),filtersize), (k==1 ? 1 : filternum*2^(k-2)) => filternum*2^(k-1) ,relu, pad=SamePad()),
    #         #x->Flux.normalise(x,dims=2),
    #         BatchNorm(filternum*2^(k-1))
    #         #x->Flux.normalise(x,dims=3)
    #         #Flux.AdaptiveMaxPool((embeddingdims,samplesize÷(poolsize*k)))
    #     ])
    # end

    append!(layers,[
        Flux.GlobalMeanPool(), 
        Flux.flatten,
        #LayerNorm(filternum*2^(numconvlayer-1),relu),
        Dense(filternum*2^(numconvlayer-1)=>nclasses)
    ])

    model = Chain(layers...)

    ps = Flux.params(model[2:end]) # taking out the EmbeddingLayer
    num_params(model) = sum(length, Flux.params(model)) 
    
    opt = ADAM(η) 
    if λ > 0 
        opt = Flux.Optimiser(WeightDecay(λ), opt)
    end
    

    @info "Batch size = $batchsize"

    @info "Dataset : $(length(train_loader.data)) train batchs | $(length(test_loader.data)) test batchs"

    @info model

    @info "initialised model with $(num_params(model[2:end])) learnable parameters."

    @info "Start Training"
    @time report(0)
    @time report(0)
    for epoch in 1:maxepoch
        @showprogress for (x, y) in train_loader
            
            gs = Flux.gradient(ps) do
                    ŷ = model(x)
                    loss(ŷ, y)
                end

            Flux.Optimise.update!(opt, ps, gs)
        end
        
        ## Printing and logging
        epoch % infotime == 0 && report(epoch)
        if checktime > 0 && epoch % checktime == 0
            !ispath(savepath) && mkpath(savepath)
            modelpath = joinpath(savepath, savename) 
            let model = cpu(model) ## return model to cpu before serialization
                BSON.@save modelpath model epoch
            end
            @info "Model saved in \"$(modelpath)\""
        end
    end
    return model
end

loss(ŷ, y) = Flux.logitcrossentropy(ŷ, y)
round4(x) = round(x, digits=4)
function eval_loss_accuracy(loader, model)
    l = 0f0
    acc = 0
    ntot = 0
    for (x, y) in loader
        
        ŷ = model(x)
        l += loss(ŷ, y) * size(x)[end]        
        acc += sum(Flux.onecold(ŷ) .== Flux.onecold(y))
        ntot += size(x)[end]
    end
    return (loss = l/ntot |> round4, acc = acc/ntot*100 |> round4)
end

function measure_accuracy(net,Xtest,Ytest,batchsize=128)
    ytest = Flux.onehotbatch(Ytest,unique(Ytest))
    Xtest=permutedims(Xtest)
    Xtest = reshape(Xtest, size(Xtest)[1], 1, size(Xtest)[2])
    test_loader = Flux.Data.DataLoader((Xtest, ytest), batchsize=batchsize, shuffle=true)
    accuracy = eval_loss_accuracy(test_loader,net)[:acc]
    println("accuracy = $accuracy")
    return accuracy
end

function plot_accuracy_v_offset(net,Xtest,Ytest,testoffsets,numbins)
    ytest = Flux.onehotbatch(Ytest,unique(Ytest))
    Xtest=permutedims(Xtest)
    Xtest = reshape(Xtest, size(Xtest)[1], 1, size(Xtest)[2])
    binoff = maximum(testoffsets)÷numbins
    bins = [binoff*i for i in 1:numbins+1]
    x=[]
    x=[[] for i in 1:length(bins)]
    @showprogress for i in eachindex(Ytest)
        for binind in eachindex(bins)
            if testoffsets[i] <= bins[binind] && testoffsets[i] > bins[binind] - binoff 
                ŷ=net(Xtest[:,:,i:i])
                
                append!(x[binind],sum(Flux.onecold(ŷ) .== Flux.onecold(ytest[:,i])))
            end
        end
    end
    Y=[]
    for l in x 
        append!(Y, sum(l)/length(l))
    end
    
    X=[binoff*i-binoff÷2 for i in 1:numbins+1]


    
    plot(X,Y;xlabel="offset from the reference sequence",ylabel="accuracy")
end
function measure_bayes_optimal_accuracy(Xtrain,Ytrain)
end
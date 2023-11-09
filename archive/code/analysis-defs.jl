using DataFrames, CSV, RCall, Statistics, MultivariateStats, Gadfly, Colors, Cairo, Fontconfig

# to use rcall
# press '$' then delete '$'
# then press '$' again
# will open R repl
# backspace to exit R repl
# more at https://juliainterop.github.io/RCall.jl/v0.12/gettingstarted.html

# note: in order to avoid a libcurl related issue
#       that occurred while importing mia i renamed
#       libcurl.so in lib/julia fldr to libcurl.so.bckp
#       this may cause issues elsewhere...

@rimport vegan

function vec2mat(dst_vec,n_sites)

    dmat = zeros(n_sites,n_sites)

    k = 0
    for i in 1:n_sites
        for j in 1:(i-1)
            k += 1        
            dmat[i,j] = dst_vec[k]
            dmat[j,i] = dst_vec[k]
        end
    end

    return dmat

end

function dst_mat(dat,mthd)

    dst = vegan.vegdist(dat,var"method"=mthd);

    dst_vec = rcopy(dst)

    n_sites = size(dat)[1]

    mat = vec2mat(dst_vec,n_sites)

    return mat

end

function mds(dat::DataFrame,mthd::String,dim::Int)

    dst = dst_mat(dat,mthd)

    M = fit(MDS,dst; distances=true, maxoutdim=dim)

    Y = predict(M)

    return Y'
    
end

# pools hosts at each time pt within a sample
function host_pool(dat,met)

    n = maximum(Int,met.Sample)

    # environment indices
    einds = met.Loc .== 0

    # host indices
    hinds = met.Loc .> 0

    # system ages recorded
    ages = unique(met.Age)

    # pooled data (empty for now)
    pdat = DataFrame()
    pmet = DataFrame()

    # pool hosts at each time pt within a sample
    for s in 1:n

        # sample indices
        sinds = met.Sample .== s

        # index of env at age a for sample s
        esinds = einds .* sinds

        # append env
        append!(pdat, dat[esinds,:])
        append!(pmet,met[esinds,["Loc","HGen","Age","Sample"]])
        
        
        # pool hosts at each time pt
        for a in ages
            
            # age indices
            ainds = met.Age .== a

            # hosts of given age indices
            inds = hinds .* ainds .* sinds

            # append pooled average rel abunds
            hpool = eachcol(dat[inds,:])
            push!(pdat,mean.(hpool))

            # append metadata
            g = mean(met[inds,:].HGen)
            push!(pmet,[1,g,a,s])

        end

    end

	return deepcopy(pdat), deepcopy(pmet)

end
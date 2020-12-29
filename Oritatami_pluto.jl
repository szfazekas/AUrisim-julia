### A Pluto.jl notebook ###
# v0.12.12

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ 1a12ca10-4391-11eb-38f3-475632625048
begin
	#using Pkg
	#Pkg.add("Colors")
	#Pkg.add("WGLMakie")
	#Pkg.add("Colors")
	using Pluto
	using PlutoUI
	using Plots
	
	using LightGraphs
	using GraphPlot

#	gr()
	using Colors
	using ColorSchemes
	using GeometryTypes
	#using GLMakie
	#using GLMakie.AbstractPlotting
	#using WGLMakie
	#using WGLMakie.AbstractPlotting
	
	#AbstractPlotting.inline!(true)
		
end

# ╔═╡ 93166c80-4861-11eb-2e9d-25075cdc8bf1
md"""
Initialization
"""


# ╔═╡ 4659a0d0-4391-11eb-08bc-5f0fa8380c40
begin
	
	#delta = 2
	#arity = 4
	
	#neighborhood = [(1,0), (1,1), (0,1), (-1,0), (-1,-1), (0,-1)]
	
	neighborhood = [[1,0], [1,1], [0,1], [-1,0], [-1,-1], [0,-1]]
	
	drawCount = 50
	
#	sigma = [1,2,1,2,3,1,2,1,3,3,2,1,1,2]
	
	#transcript = ["0","6","8","5","7","2","4","1","4","2","7","5","8","5","7","2","7","0","0"]
	
	#beads = Dict([((0,0),["0",[]]), ((-1,0),["2",[]]), ((0,1),["0",[]])])
	
	#beads = Dict([((0,0),["0",[]]), ((-1,0),["2",[]]), ((1,1),["2",[]]), ((0,1),["2",[]]), ((1,0),["2",[]]), ((-1,-1),["2",[]])])
	#Dict([((0,1),["0",[]]),((0,2),["0",[]]),((-1,1),["0",[]]),((-1,0),["0",[]])])
	
	#beads = Dict()
	
	dirs = [[1,0], [0,1], [-1,1], [-1,0], [0,-1], [1,-1]]
	
	#hood =        [          (-2,2),   (-1,2),   (0,2),
	#                     (-2,1),  (-1,1),   (0,1),    (1,1),
	#                 (-2,0),  (-1,0),   (0,0),    (1,0),   (2,0),
	#                     (-1,-1), (0,-1),   (1,-1),   (2,-1),
	#                          (0,-2),   (1,-2),   (2,-2)]
	
	
	hood =        [          [-2,2],   [-1,2],   [0,2],
	                     [-2,1],  [-1,1],   [0,1],    [1,1],
	                 [-2,0],  [-1,0],   [0,0],    [1,0],   [2,0],
	                     [-1,-1], [0,-1],   [1,-1],   [2,-1],
	                          [0,-2],   [1,-2],   [2,-2]]
	
	
	perimeter = [hood[1], hood[2], hood[3], hood[7], 
	             hood[12], hood[16], hood[19], hood[18], 
	             hood[17], hood[13], hood[8], hood[4]]
	
	initpath1 = [(0,0),(1,0), (1,1)]
	initpath2 = [(0,0),(1,0), (0,-1)]
	initpath3 = [(0,0),(1,0), (2,0), (2,1)]
	
	#currentpath = []
	#testpath = [np.array((0,1)),np.array((1,0)),np.array((1,-1)),np.array((0,-1)),np.array((-1,0)),np.array((-1,1)),np.array((0,0))]
	
	gridp = [[-2, 0], [-2, 1], [-2, 2], [-1, -1], [-1, 0], [-1, 1], [-1, 2], [0, -2], [0, -1], [0, 0], [0, 1], 
	    [0, 2], [1, -2], [1, -1], [1, 0], [1, 1], [2, -2], [2, -1], [2, 0] ]
	
	shear = [1 -0.5;0 sqrt(3)/2]
	
	#rules = Dict([("0",["0","1"]),("1",["0","1"]),("2",[])])
	
	gliderfile = "..\\AUrisim-master\\glidersample.auri.txt"
	pyramidfile = "..\\AUrisim-master\\pyramidsample.auri.txt"
	bincountfile = "..\\AUrisim-master\\binCount.auri.txt"
end

# ╔═╡ 5ef33c50-4391-11eb-0042-131969fd46d8
function loadOS(file)
    io = open(file, read=true)
    f = read(io,String)
    close(io)

	fbeads = Dict()
	
    os = split(f,"\n")

    fdelta = parse(UInt8, os[1])

    farity = parse(UInt8, os[2])

    fseed = []

    for fbeadstr in split(os[3],"->")
        fbead = split(fbeadstr,",")
        push!(fseed,[string(fbead[1]), parse(Int, fbead[2]), parse(Int, fbead[3])])
    end
    
    for fbead in fseed
        fbeads[[fbead[2],fbead[3]]] = [fbead[1],[]]
    end
    
    ftranscript = split(strip(os[4]),",")

    frules = Dict{String, Array{String}}()
    for frule in split(os[5],",")
        sides = split(frule, "=")
        if haskey(frules, sides[1])
            push!(frules[sides[1]], sides[2])
        else
            frules[sides[1]] = [sides[2]]
        end

        if haskey(frules, sides[2])
            push!(frules[sides[2]], sides[1])
        else
            frules[sides[2]] = [sides[1]]
        end
    end
    
    return [fdelta, farity, fseed, ftranscript, frules, fbeads]
end

# ╔═╡ 64b5d850-4391-11eb-0830-b72530bbf57d
function genComb(n, k)
	if n < k || k < 0
		return []
    end
	index = k
	combos = []
	comb = []

	for i=1:k
		push!(comb,i)
    end
    push!(combos,copy(comb))

	while index > 0
		if comb[index] < n-k+index
			comb[index] += 1
			for i=index+1:k
				comb[i] = comb[i-1] + 1
            end
			push!(combos,copy(comb))
			index = k
		else
			index -= 1
        end
    end
	return combos
end

# ╔═╡ 50f39900-4392-11eb-33dc-0999c97ce6b2
# generate all combinations of k elements from set, using genComb(|set|, k) above
function genCombSet(set,k)
	tmp = genComb(length(set), k)
	combos = []
	comb = []
	for i in tmp
		comb = []
		for j in i
			push!(comb,set[j])
        end
        push!(combos,copy(comb))
    end
	return combos
end

# ╔═╡ 50f5bbe0-4392-11eb-00c0-7191e55c2d31
# generate Cartesian power k of set, recursively
function genCartPower(set, k)
	cart = []
	if set == [] || k == 0
		return [[]]
    end
	tmp = genCartPower(set, k-1)
	for i in set
		for j in tmp
			push!(cart,push!(copy(j),i))
        end
    end
	return cart
end

# ╔═╡ 51076f20-4392-11eb-1ce1-012c08accf88
# generate Cartesian product of sets, recursively
function genCart(sets)
	cart = []
	if sets == [] #|| getkey(sets,[],'x')=='x'
		return [[]]
    end
	tmp = pop!(sets)
	tmp2 = genCart(sets)

	for i in tmp
		for j in tmp2
			push!(cart,append!([i],j))
        end
    end
	return cart
end

# ╔═╡ e75a6ef0-484c-11eb-3be9-7bc26e01b28e
cutoff = 5

# ╔═╡ ea695dc0-492a-11eb-3e0e-41dea86a9123
md"""
Length of transcript to try to fold: $(@bind folding_length NumberField(0:5000, default=10))
"""

# ╔═╡ f3b25b10-485e-11eb-0672-a56f14d5eaf9
md"""
Size of the beads on the plot: $(@bind beadsize Slider(1:30,default=4, show_value=true))
_______Width of plot: $(@bind canvaswidth Slider(200:100:1500, default=800, show_value=true))
_______Height of plot: $(@bind canvasheight Slider(200:100:1200, default=400, show_value=true))
"""

# ╔═╡ 583db4d2-4391-11eb-2ae9-c73fd33c55a9
function plotPath(beads, fpath, fcolored, labels=[], anim=false, fbonds=[], limit=1)
    tmp1 = []
	bonds = deepcopy(fbonds)
    minx = maxx = fpath[1][1]
    miny = maxy = fpath[1][2]
    fbeadtypes = Dict()
    fbeadtypecount = 0
    fcolorpath = []
    if fcolored == true
        fcolors = colorschemes[:tab20]
        for i=1:limit
            if fpath[i][1] < minx
                minx = fpath[i][1]
            elseif fpath[i][1]>maxx
                maxx = fpath[i][1]
            end
            if fpath[i][2] < miny
                miny = fpath[i][2]
            elseif fpath[i][2]>maxy
                maxy = fpath[i][2]
            end
        end
        if labels != []
            for i=1:limit
                if !(labels[i] in keys(fbeadtypes))
                    fbeadtypes[labels[i]] = fbeadtypecount+1
                    fbeadtypecount += 1
                end
                push!(fcolorpath, fcolors[1+mod(fbeadtypes[labels[i]],20)])
            end
        end
    end

    for i=1:limit
        push!(tmp1,shear*[fpath[i][1],fpath[i][2]])
    end
    
    tmp = [GeometryTypes.Point2f0(bead[1], bead[2]) for bead in tmp1]
	
    plt = plot(tmp, color = RGBA(0.1,0.1,0.1,0.8), linewidth = 3)
    scatter!(tmp, 
        color = fcolorpath,
        size = (canvaswidth, canvasheight),
        aspect_ratio=:equal, 
#        grid=false, 
        marker = :hexagon,
        markersize = beadsize,
        xlims = (minx-1,maxx+3),
        ylims = (miny-1,maxy+3),
#        xticks = 0:1:10,
        linealpha = 0.5,
        linewidth = 3,
        linecolor = RGBA{Float32}(0.10,0.10,0.10,1)
        )
#        linecolor = :black)
    
#	for bead in keys(beads)
#        for bond in beads[bead][2]
#            push!(bonds, [bead, bond])
#        end
#    end
    
	for bead in fpath[1:limit]
        for bond in beads[bead][2]
			if bond in fpath[1:limit]
            	push!(bonds, [bead, bond])
			end
        end
    end
	
	
    for bond in bonds
        from = GeometryTypes.Point2f0(shear*bond[1])
        to = GeometryTypes.Point2f0(shear*bond[2])

        plot!([from+(to-from)/3, to - (to-from)/3],
            color = :red,
			linewidth = 1,
			legend = :none
        )
    end
    if anim == false
        return plt
		#return bonds
    end
end

# ╔═╡ 63c99630-48b8-11eb-3750-d3afc0eb2d36
md"""
Number of beads to plot: 
$(@bind plotlimit Slider(1:folding_length; default=1, show_value=true))
"""

# ╔═╡ 158b3b60-484d-11eb-2e0d-e73b8c2fc0cd
os = loadOS(bincountfile)
#os = loadOS(pyramidfile)
#os = loadOS(gliderfile)

# ╔═╡ 321fab30-484d-11eb-3a48-bb6a1827c51d
delta = os[1]

# ╔═╡ 51179bbe-4392-11eb-3f56-9be32cb5d542
function generateDeltaPathFast(beads, path, trans)
	dpath = [[copy(path)]]
	for i=1 : delta - length(path)+1
		push!(dpath,[])
		for j in dpath[i]
			for dir in neighborhood
				if !(haskey(beads,last(j).+dir) || last(j).+dir in j)
					push!(dpath[i+1], push!(copy(j),last(j).+dir))
                end
            end
        end
    end
	return last(dpath)
end

# ╔═╡ 350984b0-484d-11eb-1d69-8725abe363a7
arity = os[2]

# ╔═╡ 384771a0-484d-11eb-1407-07c7ab19508d
seed = os[3]

# ╔═╡ 3b13d810-484d-11eb-0a83-b1982d8eb146
transcript = os[4]

# ╔═╡ 450d50e0-4856-11eb-0275-13db9e520e93
period = length(transcript)

# ╔═╡ 3d96d96e-484d-11eb-0293-b1e8656110fa
rules = os[5]

# ╔═╡ 14edea70-4980-11eb-24e1-0190280a9c1f
begin
	tmpkeys = []
	for key in keys(rules)
		push!(tmpkeys,key)
	end
	G = Graph(length(tmpkeys)) # graph with 3 vertices
	for beadtype in 1:length(tmpkeys)
		for neighbor in rules[tmpkeys[beadtype]]
			add_edge!(G,beadtype,findfirst(isequal(neighbor),tmpkeys))
		end
	end
	
# make a triangle
#add_edge!(G₁, 1, 2)
#add_edge!(G₁, 1, 3)
#add_edge!(G₁, 2, 3)

gplot(G,  layout=spectral_layout, nodelabel=tmpkeys)
end

# ╔═╡ 6bcedc90-4391-11eb-0a02-f3a177f21c2c
function valid(beads, path, bondset, sol, index, trans)
    if sol==length(bondset)*ones(length(sol))
        #print("true", path, " BONDS ", sol, "\n")
        return true
    end
	tmpBeads = Dict()
    #println(sol, index)
	for bead in path
		for dir in neighborhood
			if bead.+dir in keys(beads)
				tmpBeads[bead.+dir] = [beads[bead.+dir][1], copy(beads[bead.+dir][2])]
            end
        end
    end

	for i=2:index+1
		if !(haskey(rules,trans[i]))
			rules[trans[i]] = []
        end
		tmpBeads[path[i]] = [trans[i],[]]
		for bond in bondset[sol[i-1]]
			push!(tmpBeads[path[i]][2],path[i].+bond)
			if (haskey(tmpBeads,path[i].+bond)) && !(path[i] in tmpBeads[path[i].+bond][2])
				push!(tmpBeads[path[i].+bond][2],path[i])
				if length(tmpBeads[path[i].+bond][2]) > arity
					return false
                end
            end
        end
    end

	for i=2 : index+1
		for bond in bondset[sol[i-1]]
			if i>1 && path[i].+bond == path[i-1]
				return false
            elseif !(haskey(tmpBeads,path[i].+bond))
				return false
            elseif !(tmpBeads[path[i].+bond][1] in rules[trans[i]])
				return false
            end
        end
    end
	return true
end

# ╔═╡ 71271e00-4391-11eb-0c2e-f1fd75922269
function backtrackFastPrime(beads, path, trans, bondset)
	tmp = Dict()
	index = 1
	solutions = []
	sol = []
	bondNo = length(bondset)
	maxstrength = -1
	det = true

	for i=1:delta
		push!(sol,0)
    end

	while index > 0
		if sol[index] <= bondNo - 1
			sol[index] += 1
			strength = 0
			for i = 1:index
				strength += length(bondset[sol[i]])
            end
			if strength + (delta - index) * arity >= maxstrength
				if valid(beads, path, bondset, sol, index, trans)
					if index == delta

						if strength > maxstrength
							maxstrength = strength
							pushfirst!(solutions,deepcopy(sol))
							det = true
                        elseif strength == maxstrength
							if sol[1] == solutions[1][1]
								pushfirst!(solutions, deepcopy(sol))
							else
								det = false
                            end
                        elseif sol[1] == solutions[1][1]
							pushfirst!(solutions,deepcopy(sol))
                        end
					else
						index += 1
						sol[index] = 0
                    end
                end
            end
		else
			index -= 1
        end
    end
	return [maxstrength, deepcopy(solutions), det==true]
end

# ╔═╡ 76874eb0-4391-11eb-21bb-e52da6d0c4fd
function backtrackArity5(beads, path, trans)

    solution = []
    strength = 0
    
    
    for i=1:length(path)-1
        push!(solution,[])
    end
    
    tmpBeads = Dict()
    
    for i=2:length(path)
        tmpBeads[path[i]] = [trans[i], []]
        for dir in neighborhood
            if (path[i].+dir in keys(beads)) && (trans[i] in rules[beads[path[i].+dir][1]])
                push!(solution[i-1], dir)
                strength += 1
              
            elseif (path[i].+dir in keys(tmpBeads)) &&  (path[i].+dir != path[i-1])  && (trans[i] in rules[tmpBeads[path[i].+dir][1]])
                push!(solution[i-1], dir)
                ind = findfirst(x -> x == path[i].+dir,path)
                strength += 1
            end
        end
        
    end
    solution = [sort(sol) for sol in solution]
	return [strength, [solution], true]
end


# ╔═╡ 7bfd9f70-4391-11eb-02d7-5f825523ddec
function findFirstFast(beads, pos, trans, bondset)
	#global cutoff
	det = true
#	bondset = [[]]
    tmp = []
	maxstrength = -1
	solutions = []

    paths = generateDeltaPathFast(beads, [pos], trans)

    
    if length(paths)>100000
		print(length(paths))
		return []
    end

    if arity < cutoff
#        for i = 1 : arity
#            bondset = append!(genCombSet(neighborhood, i), bondset)
#        end
        for path in paths
            tmp = backtrackFastPrime(beads, path, trans, bondset)
            if tmp[1] > maxstrength
                if tmp[3]
                    maxstrength = tmp[1]
                    pushfirst!(solutions,[path, tmp[2], tmp[1]])
                    det = true
                else
                    det = false
                end
            elseif tmp[1] == maxstrength
                if solutions != []
                    if !(path[2] == solutions[1][1][2])
                        det = false
                    elseif !(tmp[2][1][1] == solutions[1][2][1][1])
                        det = false
                    end
                end
                pushfirst!(solutions,[path, tmp[2], tmp[1]])
            elseif tmp[1] < maxstrength
                push!(solutions,[path, tmp[2], tmp[1]])
            end
        end
        if det
            tmp = []
            for sol in solutions
                if sol[1][2] == solutions[1][1][2]
                    push!(tmp,sol)
                end
            end
            return tmp
        else
            return []
        end
    else
        for path in paths
            tmp = backtrackArity5(beads, path, trans)
            if tmp[1] > maxstrength
                maxstrength = tmp[1]
                pushfirst!(solutions,[path, tmp[2], tmp[1]])
                det = true
            elseif tmp[1] == maxstrength
                if solutions != []
                    if !(path[2] == solutions[1][1][2])
                        det = false
                    elseif !(tmp[2][1][1] == solutions[1][2][1][1])
                        det = false
                    end
                end
                pushfirst!(solutions,[path, tmp[2], tmp[1]])
                #print(solutions)
            elseif tmp[1] < maxstrength
                push!(solutions,[path, tmp[2], tmp[1]])
            end
        end
        if det
            tmp = []
            for sol in solutions
                if sol[1][2] == solutions[1][1][2]
                    push!(tmp,sol)
                end
            end
            return tmp
        else
            return []
        end
    end
end


# ╔═╡ 51206ad0-4391-11eb-3d4c-3ba8cbc17a60
function fold(fbeads, path, ftranscript, lastbead, labels, trlength)
	error = ""
	beads = deepcopy(fbeads)
	bondset = [[]]
	if arity<cutoff
		for i = 1 : arity
            bondset = append!(genCombSet(neighborhood, i), bondset)
        end
	end
		
	# FIRST STABILIZATION
    #***************************************
	#tmp = findFirstFast(path[end],append!([os[3][end][1]], ftranscript[1:delta]))
	tmp = findFirstFast(beads,path[end],append!([lastbead], ftranscript[1:delta]),bondset)
    # findFirstFast returns with an array containing elements [path,bondset-per-bead,strength-of-elongation]
    if (tmp != []) && (arity >= cutoff)
        beads[tmp[1][1][2]] = [ftranscript[1],[]]

        for bond in tmp[1][2][1][1]
                push!(beads[tmp[1][1][2]][2], tmp[1][1][2].+bond)
        end
        push!(path,tmp[1][1][2])
        push!(labels, ftranscript[1])
    
	elseif (tmp != []) && (arity<cutoff)
		beads[tmp[1][1][2]] = [ftranscript[1],[]]

        for bond in bondset[tmp[1][2][1][1]]
                push!(beads[tmp[1][1][2]][2], tmp[1][1][2].+bond)
				push!(beads[tmp[1][1][2].+bond][2], tmp[1][1][2])
        end
        push!(path,tmp[1][1][2])
        push!(labels, ftranscript[1])
	
	else
        println("Stop bead: ", 1, "  Position: ", path[end])
    end
	
	
	
	

    # STABILIZATION AFTER FIRST BEAD
    #***************************************
	
    ntranscript = cat(ftranscript,ftranscript[1:delta],dims=(1))
	
    for i=0:max(length(ntranscript)-delta, trlength)

        tmptr = ntranscript[1+mod(i,period):1+mod(i,period)+delta]

		tmp = findFirstFast(beads, path[end],tmptr,bondset)

        if (tmp != []) && (arity >= cutoff)
            beads[tmp[1][1][2]] = [tmptr[2],[]]

            for bond in tmp[1][2][1][1]
                if tmp[1][1][2].+bond in keys(beads)
                    push!(beads[tmp[1][1][2]][2], tmp[1][1][2].+bond)
					#return [bond, tmp[1]]
                end
            end
            push!(path,tmp[1][1][2])
            push!(labels, tmptr[2])
#        else
#           println("Stop bead: ", i, "  Position: ", path[end], "  Transcript window: ", tmptr)
#			break 
#        end
		
		elseif (tmp != []) && (arity<cutoff)
			beads[tmp[1][1][2]] = [tmptr[2],[]]

			for bond in bondset[tmp[1][2][1][1]]
					push!(beads[tmp[1][1][2]][2], tmp[1][1][2].+bond)
					push!(beads[tmp[1][1][2].+bond][2], tmp[1][1][2])
#					return bondset[tmp[1][2][1][1]]
			end
			push!(path,tmp[1][1][2])
			push!(labels, tmptr[2])
		else
			error = "Last stabilized bead: $(i+1)     at position: $(path[end])               Next delay-size transcript window: $(join(tmptr[2:end]))"
			break 
			#print(path[end-delta:end], "    ", transcript[i:i+delta], "\n")
			#print(beads, "\n")
		end
    end

    #plotPath2(scene, beads, path, labels)
	
	
	return [beads, path, labels,error]
end


# ╔═╡ fa566620-484d-11eb-08ea-cf27826e5d58
startbeads = deepcopy(os[6])

# ╔═╡ 46523bb0-4850-11eb-2c9b-5d7b5d828246
path = [[bead[2], bead[3]] for bead in os[3]]

# ╔═╡ 4cc0cf1e-4850-11eb-24e1-6107648b653b
labels = [bead[1] for bead in os[3]]

# ╔═╡ 8b082da0-4391-11eb-1ee3-c70c33b4cf2e
conformation = fold(startbeads, path, transcript, seed[end][1], labels, folding_length)

# ╔═╡ 860d64f0-48b4-11eb-265a-8364e73badce
plotPath(conformation[1], conformation[2], true, conformation[3], false, [], plotlimit)

# ╔═╡ 79975630-4865-11eb-1c55-290ee27b9677
html"""<style>
main {
    max-width: 1500px;
    align-self: flex-start;
    margin-left: 50px;
}
"""

# ╔═╡ 6a63bc30-4865-11eb-0396-6d4ca5f1bb07


# ╔═╡ 652808a0-4853-11eb-0dd7-cffac7eec737


# ╔═╡ 838baca0-4391-11eb-134d-e574e77f4911
function plotPath2(scene, fbeads, fpath, labels)

#    for i=1:length(fpath)
#        push!(xpath,shear*[fpath[i][1],fpath[i][2]])
#    end
    xpath = [GeometryTypes.Point2f0(shear*bead) for bead in fpath]

    fcolors = colorschemes[:tab20]
    fcolorpath = []
    fbeadtypes = Dict()
    fbeadtypecount = 0
    minx = 0
    maxx = 0
    miny = 0
    maxy = 0
    
    
    for i=1:length(fpath)
        if fpath[i][1] < minx
            minx = fpath[i][1]
        elseif fpath[i][1]>maxx
            maxx = fpath[i][1]
        end
        if fpath[i][2] < miny
            miny = fpath[i][2]
        elseif fpath[i][2]>maxy
            maxy = fpath[i][2]
        end
    end
    
    for i=1:length(xpath)
        if !(labels[i] in keys(fbeadtypes))
            fbeadtypes[labels[i]] = fbeadtypecount+1
            fbeadtypecount += 1
        end
    end


    for i=1:length(xpath)
        push!(fcolorpath, 1+mod(fbeadtypes[labels[i]],20))
    end
    
    


    lines!(scene, xpath, color = RGBA(0.1,0.1,0.1,0.8), linewidth = 3)
    scatter!(scene, xpath,
    #        Aspect(1,1),
            color = fcolors[fcolorpath],
    #        colormap = :tab20,
    #        color = :black)
    #        size = (1800, 800),
    #        aspect_ratio=:equal, 
    #        grid=false, 
            marker = :hexagon,
            markersize = 16/(abs(maxx-minx)/70))
    #        xlims = (minx-1,maxx+3),
    #        ylims = (-15,15))
    #        xticks = 0:1:10,
    #        linealpha = 0.5,
    #        linewidth = 3,
    #        linecolor = RGBA{Float32}(0.10,0.10,0.10,1)

    #ylims!(oritatami,(-10,10))
    
    bonds = []
    
    for bead in keys(fbeads)
        for bond in fbeads[bead][2]
            push!(bonds, [bead, bond])
        end
    end
    
    for bond in bonds
        from = GeometryTypes.Point2f0(shear*bond[1])
        to = GeometryTypes.Point2f0(shear*bond[2])
#        pl = plot!([tmp1[1],tmp2[1]], [tmp1[2],tmp2[2]],
        linesegments!(scene, [from+(to-from)/3, to - (to-from)/3], color = :red)
    end
#    if anim == false
#        display(scene)
#    end
    
    ylims!(scene, (miny-(maxx-minx-2*maxy+2*miny)/2, maxy + (maxx-minx-2*maxy+2*miny)/2))
    display(scene)
end

# ╔═╡ Cell order:
# ╠═1a12ca10-4391-11eb-38f3-475632625048
# ╠═14edea70-4980-11eb-24e1-0190280a9c1f
# ╟─93166c80-4861-11eb-2e9d-25075cdc8bf1
# ╟─4659a0d0-4391-11eb-08bc-5f0fa8380c40
# ╟─51206ad0-4391-11eb-3d4c-3ba8cbc17a60
# ╟─583db4d2-4391-11eb-2ae9-c73fd33c55a9
# ╟─5ef33c50-4391-11eb-0042-131969fd46d8
# ╟─64b5d850-4391-11eb-0830-b72530bbf57d
# ╟─50f39900-4392-11eb-33dc-0999c97ce6b2
# ╟─50f5bbe0-4392-11eb-00c0-7191e55c2d31
# ╟─51076f20-4392-11eb-1ce1-012c08accf88
# ╟─51179bbe-4392-11eb-3f56-9be32cb5d542
# ╟─6bcedc90-4391-11eb-0a02-f3a177f21c2c
# ╟─71271e00-4391-11eb-0c2e-f1fd75922269
# ╟─76874eb0-4391-11eb-21bb-e52da6d0c4fd
# ╟─7bfd9f70-4391-11eb-02d7-5f825523ddec
# ╟─e75a6ef0-484c-11eb-3be9-7bc26e01b28e
# ╟─321fab30-484d-11eb-3a48-bb6a1827c51d
# ╟─350984b0-484d-11eb-1d69-8725abe363a7
# ╟─384771a0-484d-11eb-1407-07c7ab19508d
# ╟─3b13d810-484d-11eb-0a83-b1982d8eb146
# ╟─3d96d96e-484d-11eb-0293-b1e8656110fa
# ╟─fa566620-484d-11eb-08ea-cf27826e5d58
# ╟─46523bb0-4850-11eb-2c9b-5d7b5d828246
# ╟─4cc0cf1e-4850-11eb-24e1-6107648b653b
# ╟─450d50e0-4856-11eb-0275-13db9e520e93
# ╟─8b082da0-4391-11eb-1ee3-c70c33b4cf2e
# ╟─ea695dc0-492a-11eb-3e0e-41dea86a9123
# ╟─f3b25b10-485e-11eb-0672-a56f14d5eaf9
# ╟─63c99630-48b8-11eb-3750-d3afc0eb2d36
# ╟─860d64f0-48b4-11eb-265a-8364e73badce
# ╠═158b3b60-484d-11eb-2e0d-e73b8c2fc0cd
# ╠═79975630-4865-11eb-1c55-290ee27b9677
# ╠═6a63bc30-4865-11eb-0396-6d4ca5f1bb07
# ╠═652808a0-4853-11eb-0dd7-cffac7eec737
# ╟─838baca0-4391-11eb-134d-e574e77f4911

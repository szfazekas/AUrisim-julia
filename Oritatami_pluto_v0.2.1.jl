### A Pluto.jl notebook ###
# v0.12.18

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
#	using Traceur
	using Pluto
	using PlutoUI
	using Plots
	using LightGraphs
	using GraphPlot
	using Colors
	using ColorSchemes
	using GeometryTypes
	md"""
	### Importing used packages
	"""
end

# ╔═╡ ea695dc0-492a-11eb-3e0e-41dea86a9123
md"""
Plotting backend: $(@bind backend Select(["gr" => "GR", "plotly" => "plotly"]))
``~~~~~~~~~~~~~~~~~~~~~~~~~``
OS definition: $(@bind osdef Select(["fromfile" => "File", "fromcell" => "Cell"]))
``~~~~~~~~~~~~~~~~~~~~~~~~~``
Length of conformation to try to fold (including the seed): $(@bind folding_length NumberField(0:5000, default=10)) beads
"""

# ╔═╡ f3b25b10-485e-11eb-0672-a56f14d5eaf9
md"""
Size of the beads on the plot: $(@bind beadsize Slider(1:30,default=4, show_value=true))
``~~~~~~~~~~~~~~~~~~`` Width of plot: $(@bind canvaswidth Slider(200:100:1500, default=800, show_value=true))
``~~~~~~~~~~~~~~~~~~`` Height of plot: $(@bind canvasheight Slider(200:100:1200, default=400, show_value=true))
"""

# ╔═╡ 63c99630-48b8-11eb-3750-d3afc0eb2d36
md"""
Number of beads to plot: 
$(@bind plotlimit Slider(1:folding_length; default=1, show_value=true))
"""

# ╔═╡ 00571140-53ee-11eb-26aa-e552dc854708
md"""
### For manual set up of OS   
Delta: $(@bind celldelta NumberField(1:12, default=3))
``~~~~~~~~~~~~~~~~~~~~~~~~~``
Arity: $(@bind cellarity NumberField(1:5, default=1))

Seed: $(@bind seedtext TextField((100,1),default="0,-1,1->0,-1,0->1,0,1->0,0,0->5,-1,-1"))

Transcript: $(@bind transcripttext TextField((100,1), default="3,0,4,1,0,5"))

Rules: $(@bind rulestext TextField((100,1),default="1=4,3=5"))
"""

# ╔═╡ 2d028b50-53f0-11eb-1f6d-9757a5b16ae4


# ╔═╡ 85206b10-5042-11eb-3cd7-0167883de0ae
md"""
### Parameters of the Oritatami system
"""

# ╔═╡ 4659a0d0-4391-11eb-08bc-5f0fa8380c40
begin
	
	#delta = 2
	#arity = 4
	#transcript = ["0","6","8","5","7","2","4","1","4","2","7","5","8","5","7","2","7","0","0"]
	#rules = Dict([("0",["0","1"]),("1",["0","1"]),("2",[])])
	
	cutoff = 5
	
	#neighborhood = [[1,0], [1,1], [0,1], [-1,0], [-1,-1], [0,-1]]
	
	#dirs = Array{Array{Int16,1}}([[1,0], [0,1], [-1,1], [-1,0], [0,-1], [1,-1]])
	#dirs = [[1,0], [0,1], [-1,1], [-1,0], [0,-1], [1,-1]]
	dirs = Vector{Vector{Int16}}([[1,0],[1,1],[0,1],[-1,0],[-1,-1],[0,-1]])
	neighborhood = dirs
	
	shear = [1 -0.5;0 sqrt(3)/2]
	
	
	gliderfile = "sampleOS\\glider.auri.txt"
	pyramidfile = "sampleOS\\pyramid.auri.txt"
	bincountfile = "sampleOS\\bincount.auri.txt"
	
	if backend == "gr"
		gr()
	else
		plotly()
	end
	
	md"""
	### Initialization of global variables
	"""
end

# ╔═╡ 65f4fca0-49b3-11eb-2ee4-072e50d481b5
	md"""
	The cell below sets up the page. If you want narrower/wider cells, change the value of max-width.
	"""

# ╔═╡ 79975630-4865-11eb-1c55-290ee27b9677
html"""<style>
main {
	max-width: 1200px;
	align-self: flex-start;
	margin-left: 50px;
}
"""

# ╔═╡ 64b5d850-4391-11eb-0830-b72530bbf57d
function genComb(n, k)
	if n < k || k <= 0
		return [[]]
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

# ╔═╡ 577d6900-4e7c-11eb-32bb-6b5cf76ccd14
md"""
### Deprecated functions, kept as fallback, if things break
- generateDeltaPathFast (replaced by __genDeltaNonRec__;  later possibly by generateDeltaTree...)
- valid (replaced by __validx__)
- backtrackFastPrime (replaced by __backtrackx__)
- backtrackArity5 (replaced by __backtrackArity5x__)
- findFirstFast (replaced by __findFirstx__  and __findNextx__)
- fold (replaced by __foldx__)
"""

# ╔═╡ 6a63bc30-4865-11eb-0396-6d4ca5f1bb07
md"""
### New data structures used
"""

# ╔═╡ c3687152-4c3d-11eb-2eab-03cc2d6cde06
mutable struct tree6
	pos :: Vector{Int16}
	btype :: String
	depth :: Int16
	children :: Vector{tree6}
	tree6() = new([Int16(0),Int16(0)], "", Int16(0), Vector{tree6}())
	tree6(pos) = new(pos, "", Int16(0), Vector{tree6}())
	tree6(pos,dp) = new(pos, "", Int16(dp), Vector{tree6}())
end

# ╔═╡ 70301ca0-4ccd-11eb-1a1b-4fcff1fc8e2c
struct stabilized
	btype :: String
	bonds :: Vector{Vector{Int16}}
	stabilized(name) = new(name,[])
	stabilized(name,bonds) = new(name,bonds)
end

# ╔═╡ d3ebb0a0-53ef-11eb-14ac-cf6c636d6dad
begin
	cellseed = []
	for fbeadstr in split(seedtext,"->")
    	fbead = split(fbeadstr,",")
        push!(cellseed,[string(fbead[1]), parse(Int16, fbead[2]), parse(Int16, fbead[3])])
	end
	
	cellbeads = Dict{Array{Int16,1},stabilized}()
	for fbead in cellseed
        cellbeads[[fbead[2],fbead[3]]] = stabilized(fbead[1])
    end
	
	celltranscript = split(strip(transcripttext),",")
	
	cellrules = Dict{String, Array{String}}()
    for frule in split(rulestext,",")
        sides = split(frule, "=")
        if haskey(cellrules, sides[1])
            push!(cellrules[sides[1]], sides[2])
        else
            cellrules[sides[1]] = [sides[2]]
        end

        if haskey(cellrules, sides[2])
            push!(cellrules[sides[2]], sides[1])
        else
            cellrules[sides[2]] = [sides[1]]
        end
    end
	md"""Setting up variables from cell ..."""
#	cellseed,celltranscript,cellrules
end

# ╔═╡ 583db4d2-4391-11eb-2ae9-c73fd33c55a9
function plotPath(beads::Dict{Array{Int16,1},stabilized}, fpath::Vector{Vector{Int16}}, fcolored::Bool, labels::Vector{String}=[], anim=false, fbonds=[], limit=1)
    tmp1 = Array{Float16,2}(undef,limit,2)
	#tmp1 = []
	tmp2 = Array{Float16,1}(undef, 2)
	from = Array{Float16,1}(undef,2)
	to = Array{Float16,1}(undef,2)
	
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
        #push!(tmp1,shear*[fpath[i][1],fpath[i][2]])
		tmp2 = shear*[fpath[i][1],fpath[i][2]]
		tmp1[i, 1] = tmp2[1]
		tmp1[i, 2] = tmp2[2]
    end
    
	plt = plot(tmp1[1:end,1], tmp1[1:end,2], color = RGBA(0.1,0.1,0.1,0.8), linewidth = 3, size = (canvaswidth, canvasheight))
    scatter!(tmp1[1:end,1], tmp1[1:end,2], 
        color = fcolorpath,
#        size = (canvaswidth, canvasheight),
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
    
	for bead in fpath[1:limit]
        for bond in beads[bead].bonds
			if bond in fpath[1:limit]
            	push!(bonds, [bead, bond])
			end
        end
    end
	
	
    for bond in bonds
        from = shear*bond[1]
        to = shear*bond[2]
        plot!([(from+(to-from)/3)[1], (to - (to-from)/3)[1]],[(from+(to-from)/3)[2], (to - (to-from)/3)[2]],
            color = :red,
			linewidth = 1,
			legend = :none
        )
    end
    if anim == false
        return plt
    end
end

# ╔═╡ 5ef33c50-4391-11eb-0042-131969fd46d8
function loadOS(file)
    io = open(file, read=true)
    f = read(io,String)
    close(io)

	fbeads = Dict{Array{Int16,1},stabilized}()
	
    os = split(f,"\n")

    fdelta = parse(UInt8, os[1])

    farity = parse(UInt8, os[2])

    fseed = []

    for fbeadstr in split(os[3],"->")
        fbead = split(fbeadstr,",")
        push!(fseed,[string(fbead[1]), parse(Int16, fbead[2]), parse(Int16, fbead[3])])
    end
    
    for fbead in fseed
        fbeads[[fbead[2],fbead[3]]] = stabilized(fbead[1])
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

# ╔═╡ 158b3b60-484d-11eb-2e0d-e73b8c2fc0cd
begin
	if osdef == "fromfile"
		#os = loadOS(bincountfile)
		#os = loadOS(pyramidfile)
		os = loadOS(gliderfile)
	else
		os = [Int16(celldelta),Int16(cellarity),cellseed,celltranscript,cellrules,cellbeads]
	end
	
	md"""
	### This cell loads the Oritatami definition file. Uncomment the one you want to use.   
	### Rerun this cell when changing the system or changing a parameter.
	"""
end

# ╔═╡ 321fab30-484d-11eb-3a48-bb6a1827c51d
delta = os[1]
#delta = 4

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
#arity = 3

# ╔═╡ 384771a0-484d-11eb-1407-07c7ab19508d
seed = os[3]

# ╔═╡ 3b13d810-484d-11eb-0a83-b1982d8eb146
transcript = Vector{String}(os[4])

# ╔═╡ 3d96d96e-484d-11eb-0293-b1e8656110fa
rules = Dict{String,Vector{String}}(os[5])

# ╔═╡ fa566620-484d-11eb-08ea-cf27826e5d58
begin
	startbeads = deepcopy(os[6])
	path = [[bead[2], bead[3]] for bead in os[3]]
	labels = [bead[1] for bead in os[3]]
	period = length(transcript)
	
	md"""
	Some internal variables store the seed in a dictionary, the path of the conformation and the labels of the vertices of the path separately for faster access
	"""
end

# ╔═╡ 6bcedc90-4391-11eb-0a02-f3a177f21c2c
function valid(beads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, bondset, sol, index::Int16, trans::Vector{String})
    if sol==length(bondset)*ones(length(sol))
        #print("true", path, " BONDS ", sol, "\n")
        return true
    end
	tmpBeads = Dict{Array{Int16,1},stabilized}()
    #println(sol, index)
	for bead in path
		for dir in neighborhood
			if bead.+dir in keys(beads)
				tmpBeads[bead.+dir] = stabilized(beads[bead.+dir].btype, copy(beads[bead.+dir].bonds))
            end
        end
    end

	for i=2:index+1
		if !(haskey(rules,trans[i]))
			rules[trans[i]] = []
        end
		tmpBeads[path[i]] = stabilized(trans[i])
		for bond in bondset[sol[i-1]]
			push!(tmpBeads[path[i]].bonds,path[i].+bond)
			if (haskey(tmpBeads,path[i].+bond)) && !(path[i] in tmpBeads[path[i].+bond].bonds)
				push!(tmpBeads[path[i].+bond].bonds,path[i])
				if length(tmpBeads[path[i].+bond].bonds) > arity
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
            elseif !(tmpBeads[path[i].+bond].btype in rules[trans[i]])
				return false
            end
        end
    end
	return true
end

# ╔═╡ 71271e00-4391-11eb-0c2e-f1fd75922269
function backtrackFastPrime(beads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, trans::Vector{String}, bondset)
#	tmp = Dict()
	index = Int16(1)
	solutions = []
	sol = []
	bondNo = length(bondset)
	maxstrength = Int16(-1)
	det = true

	for i=1:delta
		push!(sol,Int16(0))
    end

	while index > 0
		if sol[index] <= bondNo - 1
			sol[index] += 1
			strength = Int16(0)
			for i = 1:index
				strength += length(bondset[sol[i]])
            end
			if strength + (delta - index) * arity >= maxstrength
				if valid(beads, path, bondset, sol, Int16(index), trans)
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
						sol[index] = Int16(0)
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
function backtrackArity5(beads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, trans::Vector{String})

    solution = []
    strength = Int16(0)
    
    
    for i=1:length(path)-1
        push!(solution,[])
    end
    
    tmpBeads = Dict{Array{Int16,1},stabilized}()
    
    for i=2:length(path)
        tmpBeads[path[i]] = stabilized(trans[i])
        for dir in neighborhood
            if (path[i].+dir in keys(beads)) && (trans[i] in rules[beads[path[i].+dir].btype])
                push!(solution[i-1], dir)
                strength += 1
              
            elseif (path[i].+dir in keys(tmpBeads)) &&  (path[i].+dir != path[i-1])  && (trans[i] in rules[tmpBeads[path[i].+dir].btype])
                push!(solution[i-1], dir)
                ind = findfirst(x -> x == path[i].+dir,path)
                strength += 1
            end
        end
        
    end
    solution = [sort(sol) for sol in solution]
	return [strength, [solution], true]
end


# ╔═╡ 245fcfc0-4dcf-11eb-07b3-45899458d079
mutable struct bondcombo
	bonds :: Vector{Vector{Vector{Vector{Int16}}}}
	strength :: Int16
	det :: Bool
#	bondcombo() = new(Vector{Vector{Int16}}(undef,delta),0, true)
end

# ╔═╡ 90a0b682-4dd0-11eb-1e1b-af59989495e9
struct elongation
	path :: Vector{Vector{Int16}}
	bondseq :: bondcombo
end

# ╔═╡ e5d4d620-4e3b-11eb-1168-098a4f0be84d
struct elongtocheck
	el :: elongation
	neighbors :: Vector{Vector{Int16}}
	elongtocheck(elpar) = new(elpar,Vector{Vector{Int16}}())
end

# ╔═╡ 5c070640-4e70-11eb-2fb9-192fb5be5aee
md"""
### The cleaner code for folding. Not yet faster than the original, but perhaps easier to follow/improve
"""

# ╔═╡ 2b8bb7e0-4c9f-11eb-1e17-c173091e604a
function genDeltaNonRec(beads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, trans::Array{String,1})
	paths = Vector{Vector{Vector{Int16}}}()
	pos = zeros(Int16,2)
	tmp = zeros(Int16,2)
	current = zeros(Int16,delta+length(path))
	currentpath = Vector{Vector{Int16}}(undef,delta+length(path)) #where {T::Vector{Int16}}
	fill!(currentpath,9999*ones(Int16,2))
	currentpath[1:length(path)]=path
	index = length(path)+1
	current[index] = 0
	pos = path[1]
	pathlength=length(path)
	lengthbound = delta + pathlength
	
	while index > pathlength
		if current[index] < 6
			current[index] += 1
			if index > pathlength
				tmp = currentpath[index-1].+dirs[current[index]]
			else
				tmp = pos.+dirs[current[index]]
			end
			#pos = pos.+dirs[current[index]] 
			if !(tmp in keys(beads)) && !(tmp in currentpath)
				currentpath[index] = copy(tmp)
				if index == lengthbound
					push!(paths, copy(currentpath))
				else
					index += 1
					current[index] = Int16(0)
				end
			end
		else
			currentpath[index] = 9999*ones(Int16,2)
			index -= 1
		end
	end
	return paths
end

# ╔═╡ 7bfd9f70-4391-11eb-02d7-5f825523ddec
function findFirstFast(beads::Dict{Array{Int16,1},stabilized}, pos::Vector{Int16}, trans::Vector{String}, bondset)
	#global cutoff
	det = true
#	bondset = [[]]
    tmp = []
	maxstrength = Int16(-1)
	solutions = []

    #paths = generateDeltaPathFast(beads, [pos], trans)
	paths = genDeltaNonRec(beads, [pos], trans)
    
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
function fold(fbeads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, ftranscript::Vector{String}, lastbead::String, labels::Vector{String}, trlength::Int64)
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

	tmp = findFirstFast(beads,path[end],append!([lastbead], ftranscript[1:delta]),bondset)
    # findFirstFast returns with an array containing elements [path,bondset-per-bead,strength-of-elongation]

	if (tmp != []) && (arity >= cutoff)
        beads[tmp[1][1][2]] = stabilized(ftranscript[1])

        for bond in tmp[1][2][1][1]
                push!(beads[tmp[1][1][2]].bonds, tmp[1][1][2].+bond)
        end
        push!(path,tmp[1][1][2])
        push!(labels, ftranscript[1])
    
	elseif (tmp != []) && (arity<cutoff)
		beads[tmp[1][1][2]] = stabilized(ftranscript[1])

        for bond in bondset[tmp[1][2][1][1]]
                push!(beads[tmp[1][1][2]].bonds, tmp[1][1][2].+bond)
				push!(beads[tmp[1][1][2].+bond].bonds, tmp[1][1][2])
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
            beads[tmp[1][1][2]] = stabilized(tmptr[2])

            for bond in tmp[1][2][1][1]
                if tmp[1][1][2].+bond in keys(beads)
                    push!(beads[tmp[1][1][2]].bonds, tmp[1][1][2].+bond)
                end
            end
            push!(path,tmp[1][1][2])
            push!(labels, tmptr[2])
		
		elseif (tmp != []) && (arity<cutoff)
			beads[tmp[1][1][2]] = stabilized(tmptr[2])

			for bond in bondset[tmp[1][2][1][1]]
					push!(beads[tmp[1][1][2]].bonds, tmp[1][1][2].+bond)
					push!(beads[tmp[1][1][2].+bond].bonds, tmp[1][1][2])
			end
			push!(path,tmp[1][1][2])
			push!(labels, tmptr[2])

		else
			error = "Last stabilized bead: $(i+1)     at position: $(path[end])               Next delay-size transcript window: $(join(tmptr[2:end]))"
			break 
		end
    end
	
	return [beads, path, labels,error]
end


# ╔═╡ 2e35b6d0-4ddf-11eb-1ce7-b5f062c3c4eb
function validx(beads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, solspace::Vector{Vector{Vector{Vector{Int16}}}}, sol::Vector{Int16}, index::Int16, trans::Vector{String})::Bool
	tmpbeads = deepcopy(beads)
	for i in 2:1+index
#		if length(beads[bead].bonds) > arity
#			return false
#		end
		tmpbeads[path[i]]=stabilized(trans[i])
		for neighbor in solspace[i-1][sol[i-1]]
			push!(tmpbeads[path[i]].bonds,neighbor)
			push!(tmpbeads[neighbor].bonds, path[i])
			if length(tmpbeads[neighbor].bonds) > arity
				return false
			end
		end
#		for neighbor in tmpbeads[path[i]].bonds
#			if length(beads[neighbor].bonds) > arity) || !(path[i] in beads[neighbor].bonds)
#				return false
#			end
#		end
	end
	return true
end

# ╔═╡ 2db404b0-4dd4-11eb-0811-4ffe0855dc4f
function backtrackx(beads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, trans::Vector{String})::bondcombo
	index = Int16(1)
	solutions = Vector{Vector{Vector{Vector{Int16}}}}()
	neighbors = Vector{Vector{Vector{Int16}}}(undef,delta)
	sol = zeros(Int16,delta)
	maxstrength = Int16(-1)
	det = true
	solspace = Vector{Vector{Vector{Vector{Int16}}}}(undef,delta)
	tmpbeads = deepcopy(beads)
	
	for i=2:delta+1
		tmpbeads[path[i]] = stabilized(trans[i])
		if !haskey(rules,trans[i])
			rules[trans[i]] = Vector{String}()
		end
	end
	
	for i=2:delta+1
		neighbors[i-1] = []
		for dir in dirs
			if haskey(tmpbeads, path[i].+dir) && 
			(!(i>1 && path[i-1]==path[i].+dir) && !(i<delta+1 && path[i+1]==path[i].+dir))
				if tmpbeads[path[i].+dir].btype in rules[trans[i]]
					push!(neighbors[i-1],path[i].+dir)
				end
			end
		end
		solspace[i-1] = []
		for j in 0:arity
			append!(solspace[i-1], genCombSet(neighbors[i-1],j))
		end
    end

	while index > 0
		if sol[index] < length(solspace[index])
			sol[index] += 1
			strength = Int16(0)
			for i = 1:index
				strength += length(solspace[i][sol[i]])
            end
			if strength + (delta - index) * arity >= maxstrength
				if validx(tmpbeads, path, solspace, sol, index, trans)
					if index == delta

						if strength > maxstrength
							maxstrength = strength
							pushfirst!(solutions, deepcopy(map(x->solspace[x][sol[x]],1:delta)))
							det = true
                        elseif strength == maxstrength
							if solspace[1][sol[1]] == solutions[1][1]
								pushfirst!(solutions, deepcopy(map(x->solspace[x][sol[x]],1:delta)))
							else
								det = false
                            end
                        elseif solspace[1][sol[1]] == solutions[1][1]
							push!(solutions,copy(map(x->solspace[x][sol[x]],1:delta)))
                        end
					else
						index += Int16(1)
						sol[index] = Int16(0)
                    end
                end
            end
		else
			index -= Int16(1)
        end
    end
	return bondcombo(solutions, maxstrength, det)
end

# ╔═╡ 4f127f70-4e73-11eb-3a5c-933eb957e069
function backtrackArity5x(beads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, trans::Vector{String})

    solution = []
    strength = Int16(0)
    
    
    for i=1:length(path)-1
        push!(solution,[])
    end
    
    tmpBeads = Dict{Array{Int16,1},stabilized}()
    
    for i=2:length(path)
        tmpBeads[path[i]] = stabilized(trans[i])
        for dir in dirs
            if (path[i].+dir in keys(beads)) && (trans[i] in rules[beads[path[i].+dir].btype])
                push!(solution[i-1], path[i].+dir)
                strength += 1
              
            elseif (path[i].+dir in keys(tmpBeads)) &&  (path[i].+dir != path[i-1])  && (trans[i] in rules[tmpBeads[path[i].+dir].btype])
                push!(solution[i-1], path[i].+dir)
                ind = findfirst(x -> x == path[i].+dir,path)
                strength += 1
            end
        end
        
    end
    solution = [sort(sol) for sol in solution]
	return bondcombo([solution], strength, true)
end

# ╔═╡ b1f70610-4dce-11eb-19bb-79506dcc9c9d
function findFirstx(beads::Dict{Array{Int16,1},stabilized}, pos::Vector{Int16}, trans::Vector{String})::Vector{elongation}
	#global cutoff
	det = true
    tmp = []
	maxstrength = Int16(-1)
	solutions = Vector{elongation}()

    #paths = generateDeltaPathFast(beads, [pos], trans)
	paths = genDeltaNonRec(beads, [pos], trans)
    
    if length(paths)>10000000
		print(length(paths))
		return []
    end

    if arity < cutoff
        for path in paths
            tmp = backtrackx(beads, path, trans)
            if tmp.strength > maxstrength
                if tmp.det
                    maxstrength = tmp.strength
                    pushfirst!(solutions,elongation(path,deepcopy(tmp)))
                    det = true
                else
                    det = false
                end
            elseif tmp.strength == maxstrength
                if solutions != []
                    if !(path[2] == solutions[1].path[2])
                        det = false
                    elseif !(tmp.bonds[1][1] == solutions[1].bondseq.bonds[1][1])
                        det = false
                    end
                end
                pushfirst!(solutions,elongation(path,deepcopy(tmp)))
            elseif tmp.strength < maxstrength
                push!(solutions,elongation(path,deepcopy(tmp)))
            end
        end
        if det
            return solutions[map(x->x.path[2]==solutions[1].path[2],solutions)]
        else
            return []#solutions
        end
    else
        for path in paths
            tmp = backtrackArity5x(beads, path, trans)
            if tmp.strength > maxstrength
                maxstrength = tmp.strength
                pushfirst!(solutions,elongation(path, deepcopy(tmp)))
                det = true
            elseif tmp.strength == maxstrength
                if solutions != []
                    if !(path[2] == solutions[1].path[2])
                        det = false
                    elseif !(tmp.bonds[1][1] == solutions[1].bondseq.bonds[1][1])
                        det = false
                    end
                end
                pushfirst!(solutions,elongation(path,deepcopy(tmp)))
            elseif tmp.strength < maxstrength
                push!(solutions,elongation(path,deepcopy(tmp)))
            end
        end
        if det
            return solutions[map(x->x.path[2]==solutions[1].path[2],solutions)]
        else
            return []#solutions
        end
    end
end


# ╔═╡ f53c6de0-4e3a-11eb-365c-591d4b585e05
function findNextx(beads::Dict{Array{Int16,1},stabilized}, elongs::Vector{elongation}, trans::Vector{String})::Vector{elongation}
	det = true
	
	maxstrength = Int16(-1)
	solutions = Vector{elongation}()
	paths = Vector{Vector{Vector{Int16}}}()
#	xpaths = Vector{elongtocheck}()

    

    if arity < cutoff
		for elong in elongs
			for dir in dirs
				if !haskey(beads, elong.path[end].+dir) && !(elong.path[end].+dir in elong.path)
					push!(paths, vcat(elong.path[2:end],[elong.path[end].+dir]))
#					push!(xpaths, elongtocheck(vcat(elong.path[2:end],[elong.path[end].+dir]),elong.bondseq)
#					for nb in dirs
#						if haskey(beads,elong.path[end].+dir.+nb) && (dir.+nb != [0,0])
#							push!(xpaths[end].neighbors, elong.path[end].+dir.+nb)
#						endd
#					end
				end
			end
		end
		
				
				
				
#		for elong in xpaths
#			tmp :: bondcombo = backtrackone(beads, elong, trans, elong.neighbors)
#			if tmp.strength > maxstrength
#				if tmp.det
#					maxstrength = tmp.strength
#					pushfirst!(solutions,elongation(copy(path),deepcopy(tmp)))
#					det = true
#				else
#					det = false
#				end
#			elseif tmp.strength == maxstrength
#				if solutions != []
#					if path[2] != solutions[1].path[2]
#						det = false
#					elseif !(tmp.bonds[1][1] == solutions[1].bondseq.bonds[1][1])
#						det = false
#					end
#				end
#				pushfirst!(solutions,elongation(copy(path),deepcopy(tmp)))
#			elseif tmp.strength < maxstrength
#				push!(solutions,elongation(copy(path),deepcopy(tmp)))
#			end
#		end
		
#		if det
#			return solutions[map(x->x.path[2]==solutions[1].path[2],solutions)]
#		else
#			return []#solutions
#		end
				
		for path in paths
			tmp :: bondcombo = backtrackx(beads, path, trans)
			if tmp.strength > maxstrength
				if tmp.det
					maxstrength = tmp.strength
					pushfirst!(solutions,elongation(copy(path),deepcopy(tmp)))
					det = true
				else
					det = false
				end
			elseif tmp.strength == maxstrength
				if solutions != []
					if path[2] != solutions[1].path[2]
						det = false
					elseif false#!(solutions[1].bondseq.bonds[1][1] == tmp.bonds[1][1]) 
						det = false
					end
				end
				pushfirst!(solutions,elongation(copy(path),deepcopy(tmp)))
			elseif tmp.strength < maxstrength
				push!(solutions,elongation(copy(path),deepcopy(tmp)))
			end
		end
		
		if det
			return solutions[map(x->x.path[2]==solutions[1].path[2],solutions)]
		else
			return []#solutions
		end
    else
        for path in paths
            tmp = backtrackArity5x(beads, path, trans)
            if tmp.strength > maxstrength
                maxstrength = tmp.strength
                pushfirst!(solutions,elongation(path, deepcopy(tmp)))
                det = true
            elseif tmp.strength == maxstrength
                if solutions != []
                    if !(path[2] == solutions[1].path[2])
                        det = false
                    elseif !(tmp.bonds[1][1] == solutions[1].bondseq.bonds[1][1])
                        det = false
                    end
                end
                pushfirst!(solutions,elongation(path,deepcopy(tmp)))
            elseif tmp.strength < maxstrength
                push!(solutions,elongation(path,deepcopy(tmp)))
            end
        end
        if det
            return solutions[map(x->x.path[2]==solutions[1].path[2],solutions)]
        else
            return solutions
        end
    end
end

# ╔═╡ f8a981c0-4de5-11eb-237c-fd23db560e94
function foldx(fbeads::Dict{Array{Int16,1},stabilized}, path::Vector{Vector{Int16}}, ftranscript::Vector{String}, lastbead::String, labels::Vector{String}, trlength::Int64)
	error = ""
	beads = deepcopy(fbeads)
#	bondset = [[]]
		
	# FIRST STABILIZATION
    #***************************************

	tmp = findFirstx(beads,path[end],append!([lastbead], ftranscript[1:delta]))
    # findFirstFast returns with an array containing elements [path,bondset-per-bead,strength-of-elongation]

#******************
#	if (tmp != []) && (arity >= cutoff)
#        beads[tmp[1][1][2]] = stabilized(ftranscript[1])

#        for bond in tmp[1][2][1][1]
#                push!(beads[tmp[1][1][2]].bonds, tmp[1][1][2].+bond)
#        end
#        push!(path,tmp[1][1][2])
#        push!(labels, ftranscript[1])
#    
#	elseif (tmp != []) && (arity<cutoff)
#******************
	if (tmp != []) #&& (arity<cutoff)
		beads[tmp[1].path[2]] = stabilized(ftranscript[1])

        for bond in tmp[1].bondseq.bonds[1][1]
			if haskey(beads, bond)
                push!(beads[tmp[1].path[2]].bonds, bond)
				push!(beads[bond].bonds, tmp[1].path[2])
			end
        end
        push!(path,tmp[1].path[2])
        push!(labels, ftranscript[1])
	
	else
        println("Stop bead: ", 1, "  Position: ", path[end])
    end
	
	
	
	

    # STABILIZATION AFTER FIRST BEAD
    #***************************************
	
    ntranscript = cat(ftranscript,ftranscript[1:delta],dims=(1))
	
    for i=0:max(length(ntranscript)-delta, trlength)

        tmptr = ntranscript[1+mod(i,period):1+mod(i,period)+delta]

		if arity < cutoff 
			tmp = findNextx(beads, deepcopy(tmp), tmptr)
		else
			tmp = findFirstx(beads, copy(tmp[1].path[2]), tmptr)
		end
#******************
#        if (tmp != []) && (arity >= cutoff)
#            beads[tmp[1][1][2]] = stabilized(tmptr[2])

#            for bond in tmp[1][2][1][1]
#                if tmp[1][1][2].+bond in keys(beads)
#                    push!(beads[tmp[1][1][2]].bonds, tmp[1][1][2].+bond)
#                end
#            end
#            push!(path,tmp[1][1][2])
#            push!(labels, tmptr[2])

#		elseif (tmp != []) && (arity<cutoff)
#*******************		
		if (tmp != []) #&& (arity<cutoff)
			beads[tmp[1].path[2]] = stabilized(tmptr[2])

			for bond in tmp[1].bondseq.bonds[1][1]
				if haskey(beads, bond)
					push!(beads[tmp[1].path[2]].bonds, bond)
					push!(beads[bond].bonds, tmp[1].path[2])
				end
			end
			push!(path,tmp[1].path[2])
			push!(labels, tmptr[2])

		else
			error = "Last stabilized bead: $(i+1)     at position: $(path[end])               Next delay-size transcript window: $(join(tmptr[2:end]))"
			break 
		end
    end
	
	return [beads, path, labels,error]
end

# ╔═╡ 8b082da0-4391-11eb-1ee3-c70c33b4cf2e
conformation = foldx(startbeads, path, transcript, seed[end][1], labels, folding_length);


# ╔═╡ 860d64f0-48b4-11eb-265a-8364e73badce
plotPath(conformation[1], conformation[2], true, conformation[3], false, [], plotlimit)

# ╔═╡ 92d8d6b0-49b3-11eb-2b82-07a4ab8dcaee
md"""
### Some experimental stuff will be in the cells below, like searching for the optimal conformation in a tree rather than array of paths, alternative plotting and visualization for the rule graph. Not yet working ...
"""

# ╔═╡ a2193400-4e6b-11eb-07f5-dd856c6f5dc7
function backtrackone(beads::Dict{Array{Int16,1},stabilized}, elong::elongation, trans::Vector{String}, neighbors::Vector{Vector{Int16}}) :: bondcombo
	maxstrength = 0
	for i in 0:arity
		for lastcombo in genCombSet(neighbors,i)
			for combo in elong.bondseq
				for bond in lastcombo
					# add path[end] to the bondlist of "bond" if possible (arity is not exceeded for "bond")
					# compute strength ...
				end
			end
		end
	end
end

# ╔═╡ 70a71be0-4e42-11eb-1cdb-e5592529a943
genDeltaNonRec(startbeads,Vector{Vector{Int16}}([[-1,-1]]), ["5","3","0","4","1"])

# ╔═╡ 57b858f0-4de9-11eb-3128-c1dded54a2e3
#tmp = findFirstx(startbeads,Vector{Int16}([-1,-1]),["5","3","0","4","1","0","5"],[])

# ╔═╡ de600e92-4e41-11eb-1290-b3c24b9bcfd7
xbeads = deepcopy(startbeads)

# ╔═╡ df8308c0-4e57-11eb-0fd3-0fd5d532a2b9
xbeads

# ╔═╡ 12c8ef80-4c5d-11eb-0e3e-9f50ee5fb253
deltap=4

# ╔═╡ 68f3c870-4c4a-11eb-19bb-2198d13c7f53
function expand(beads::Dict{Array{Int16,1},stabilized}, t::tree6, depth::Int16, path::Array{Array{Int16,1},1}, trans::Array{String,1})
	if depth == deltap
		return nothing
	end
	tmppath = copy(path)
	for dir in dirs
		if !(t.pos.+dir in keys(beads)) && !(t.pos.+dir in tmppath)
			push!(t.children, tree6(t.pos.+dir, t.depth+1))
			push!(tmppath, t.pos)
			t.children[end].btype = trans[1]
			expand(beads, t.children[end],  Int16(depth+1), tmppath, trans[2:end])
		end
	end
	return nothing
end

# ╔═╡ 57513d30-4c3d-11eb-1e04-a97936ea3884
function generateDeltaTree(beads::Dict{Array{Int16,1},stabilized}, pos::Array{Int16,1}, trans::Array{String,1})
	deltatree = tree6(pos)
	deltatree.btype = beads[pos].btype
	#tmpbeads = [""]
	for dir in dirs
		if !(deltatree.pos.+dir in keys(beads))
			push!(deltatree.children,tree6(deltatree.pos.+dir, deltatree.depth+1))
			deltatree.children[end].btype = trans[1]
			expand(beads, deltatree.children[end], deltatree.children[end].depth, [deltatree.pos], trans[2:end])
		end
	end
	return deltatree
end

# ╔═╡ 2cf7ccae-4dec-11eb-3c10-7f6022d1f716
#backtrackx(startbeads,Vector{Vector{Int16}}([[-1,-1],[0,-1],[1,0],[1,1]]),["5","3","0","4","1","0","4"],[])

# ╔═╡ 1c811760-4c57-11eb-0631-2bf35a5ac258
dpaths = generateDeltaTree(Dict([Int16(0),Int16(0)]=>stabilized("0")),[Int16(0),Int16(0)], ["B1","B2", "B3","B4","","","","","","","","","","","","","",""]);

# ╔═╡ df9785a0-4cab-11eb-3b40-03868856abc5
dpaths3= genDeltaNonRec(Dict(zeros(Int16,2)=>stabilized("0")),[zeros(Int16,2)], ["B1","B2", "B3","B4","","","","","","","","","","","","","",""]);

# ╔═╡ 2c372c60-4c9a-11eb-1103-03d32d695b90
#(Base.summarysize(dpaths), Base.summarysize(dpaths2), Base.summarysize(dpaths3))

# ╔═╡ Cell order:
# ╟─ea695dc0-492a-11eb-3e0e-41dea86a9123
# ╟─f3b25b10-485e-11eb-0672-a56f14d5eaf9
# ╟─63c99630-48b8-11eb-3750-d3afc0eb2d36
# ╠═860d64f0-48b4-11eb-265a-8364e73badce
# ╠═158b3b60-484d-11eb-2e0d-e73b8c2fc0cd
# ╟─00571140-53ee-11eb-26aa-e552dc854708
# ╟─d3ebb0a0-53ef-11eb-14ac-cf6c636d6dad
# ╠═2d028b50-53f0-11eb-1f6d-9757a5b16ae4
# ╟─85206b10-5042-11eb-3cd7-0167883de0ae
# ╟─321fab30-484d-11eb-3a48-bb6a1827c51d
# ╟─350984b0-484d-11eb-1d69-8725abe363a7
# ╟─384771a0-484d-11eb-1407-07c7ab19508d
# ╟─3b13d810-484d-11eb-0a83-b1982d8eb146
# ╟─3d96d96e-484d-11eb-0293-b1e8656110fa
# ╟─fa566620-484d-11eb-08ea-cf27826e5d58
# ╠═8b082da0-4391-11eb-1ee3-c70c33b4cf2e
# ╟─1a12ca10-4391-11eb-38f3-475632625048
# ╟─4659a0d0-4391-11eb-08bc-5f0fa8380c40
# ╟─65f4fca0-49b3-11eb-2ee4-072e50d481b5
# ╠═79975630-4865-11eb-1c55-290ee27b9677
# ╟─583db4d2-4391-11eb-2ae9-c73fd33c55a9
# ╠═5ef33c50-4391-11eb-0042-131969fd46d8
# ╟─64b5d850-4391-11eb-0830-b72530bbf57d
# ╟─50f39900-4392-11eb-33dc-0999c97ce6b2
# ╟─50f5bbe0-4392-11eb-00c0-7191e55c2d31
# ╟─51076f20-4392-11eb-1ce1-012c08accf88
# ╟─577d6900-4e7c-11eb-32bb-6b5cf76ccd14
# ╟─51179bbe-4392-11eb-3f56-9be32cb5d542
# ╟─6bcedc90-4391-11eb-0a02-f3a177f21c2c
# ╟─71271e00-4391-11eb-0c2e-f1fd75922269
# ╟─76874eb0-4391-11eb-21bb-e52da6d0c4fd
# ╟─7bfd9f70-4391-11eb-02d7-5f825523ddec
# ╟─51206ad0-4391-11eb-3d4c-3ba8cbc17a60
# ╟─6a63bc30-4865-11eb-0396-6d4ca5f1bb07
# ╠═c3687152-4c3d-11eb-2eab-03cc2d6cde06
# ╠═70301ca0-4ccd-11eb-1a1b-4fcff1fc8e2c
# ╠═245fcfc0-4dcf-11eb-07b3-45899458d079
# ╠═90a0b682-4dd0-11eb-1e1b-af59989495e9
# ╠═e5d4d620-4e3b-11eb-1168-098a4f0be84d
# ╠═5c070640-4e70-11eb-2fb9-192fb5be5aee
# ╟─2b8bb7e0-4c9f-11eb-1e17-c173091e604a
# ╟─2e35b6d0-4ddf-11eb-1ce7-b5f062c3c4eb
# ╟─2db404b0-4dd4-11eb-0811-4ffe0855dc4f
# ╟─4f127f70-4e73-11eb-3a5c-933eb957e069
# ╟─b1f70610-4dce-11eb-19bb-79506dcc9c9d
# ╟─f53c6de0-4e3a-11eb-365c-591d4b585e05
# ╟─f8a981c0-4de5-11eb-237c-fd23db560e94
# ╟─92d8d6b0-49b3-11eb-2b82-07a4ab8dcaee
# ╟─68f3c870-4c4a-11eb-19bb-2198d13c7f53
# ╟─57513d30-4c3d-11eb-1e04-a97936ea3884
# ╟─a2193400-4e6b-11eb-07f5-dd856c6f5dc7
# ╟─70a71be0-4e42-11eb-1cdb-e5592529a943
# ╠═57b858f0-4de9-11eb-3128-c1dded54a2e3
# ╟─de600e92-4e41-11eb-1290-b3c24b9bcfd7
# ╟─df8308c0-4e57-11eb-0fd3-0fd5d532a2b9
# ╟─12c8ef80-4c5d-11eb-0e3e-9f50ee5fb253
# ╠═2cf7ccae-4dec-11eb-3c10-7f6022d1f716
# ╠═1c811760-4c57-11eb-0631-2bf35a5ac258
# ╠═df9785a0-4cab-11eb-3b40-03868856abc5
# ╠═2c372c60-4c9a-11eb-1103-03d32d695b90

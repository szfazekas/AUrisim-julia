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
	using Pluto
	using PlutoUI
	using Plots
	using Colors
	using ColorSchemes
	md"""
	Importing used packages
	"""
end

# ╔═╡ ea695dc0-492a-11eb-3e0e-41dea86a9123
md"""
Length of conformation to try to fold (including the seed): $(@bind folding_length NumberField(0:5000, default=10))
_          _
Plotting backend: $(@bind backend Select(["gr" => "GR", "plotly" => "plotly"]))
"""

# ╔═╡ f3b25b10-485e-11eb-0672-a56f14d5eaf9
md"""
Size of the beads on the plot: $(@bind beadsize Slider(1:30,default=4, show_value=true))
_______Width of plot: $(@bind canvaswidth Slider(200:100:1500, default=800, show_value=true))
_______Height of plot: $(@bind canvasheight Slider(200:100:1200, default=400, show_value=true))
"""

# ╔═╡ 63c99630-48b8-11eb-3750-d3afc0eb2d36
md"""
Number of beads to plot: 
$(@bind plotlimit Slider(1:folding_length; default=1, show_value=true))
"""

# ╔═╡ 07cf29ee-4b31-11eb-3740-07c6145ffafd
if backend == "gr"
	gr()
else
	plotly()
end

# ╔═╡ 4659a0d0-4391-11eb-08bc-5f0fa8380c40
begin
	
	#delta = 2
	#arity = 4
	#transcript = ["0","6","8","5","7","2","4","1","4","2","7","5","8","5","7","2","7","0","0"]
	#rules = Dict([("0",["0","1"]),("1",["0","1"]),("2",[])])
	
	cutoff = 5
	
	neighborhood = [[1,0], [1,1], [0,1], [-1,0], [-1,-1], [0,-1]]
	
	dirs = [[1,0], [0,1], [-1,1], [-1,0], [0,-1], [1,-1]]
	
	shear = [1 -0.5;0 sqrt(3)/2]
	
	
	gliderfile = "sampleOS\\glider.auri.txt"
	pyramidfile = "sampleOS\\pyramid.auri.txt"
	bincountfile = "sampleOS\\bincount.auri.txt"
	
	md"""
	Initialization of global variables
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

# ╔═╡ 583db4d2-4391-11eb-2ae9-c73fd33c55a9
function plotPath(beads, fpath, fcolored, labels=[], anim=false, fbonds=[], limit=1)
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
        for bond in beads[bead][2]
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

# ╔═╡ 158b3b60-484d-11eb-2e0d-e73b8c2fc0cd
os = loadOS(bincountfile)
#os = loadOS(pyramidfile)
#os = loadOS(gliderfile)

# ╔═╡ 321fab30-484d-11eb-3a48-bb6a1827c51d
delta = os[1]

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

# ╔═╡ fa566620-484d-11eb-08ea-cf27826e5d58
startbeads = deepcopy(os[6])

# ╔═╡ 46523bb0-4850-11eb-2c9b-5d7b5d828246
path = [[bead[2], bead[3]] for bead in os[3]]

# ╔═╡ 4cc0cf1e-4850-11eb-24e1-6107648b653b
labels = [bead[1] for bead in os[3]]

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
                end
            end
            push!(path,tmp[1][1][2])
            push!(labels, tmptr[2])
		
		elseif (tmp != []) && (arity<cutoff)
			beads[tmp[1][1][2]] = [tmptr[2],[]]

			for bond in bondset[tmp[1][2][1][1]]
					push!(beads[tmp[1][1][2]][2], tmp[1][1][2].+bond)
					push!(beads[tmp[1][1][2].+bond][2], tmp[1][1][2])
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


# ╔═╡ 8b082da0-4391-11eb-1ee3-c70c33b4cf2e
conformation = fold(startbeads, path, transcript, seed[end][1], labels, folding_length)

# ╔═╡ 860d64f0-48b4-11eb-265a-8364e73badce
plotPath(conformation[1], conformation[2], true, conformation[3], false, [], plotlimit)

# ╔═╡ 6a63bc30-4865-11eb-0396-6d4ca5f1bb07


# ╔═╡ 92d8d6b0-49b3-11eb-2b82-07a4ab8dcaee
md"""
Some experimental stuff will be in the cells below, like an alternative plotting and visualization for the rule graph. Not yet working ...
"""

# ╔═╡ Cell order:
# ╟─ea695dc0-492a-11eb-3e0e-41dea86a9123
# ╟─f3b25b10-485e-11eb-0672-a56f14d5eaf9
# ╟─63c99630-48b8-11eb-3750-d3afc0eb2d36
# ╠═860d64f0-48b4-11eb-265a-8364e73badce
# ╠═158b3b60-484d-11eb-2e0d-e73b8c2fc0cd
# ╟─321fab30-484d-11eb-3a48-bb6a1827c51d
# ╟─350984b0-484d-11eb-1d69-8725abe363a7
# ╟─384771a0-484d-11eb-1407-07c7ab19508d
# ╟─3b13d810-484d-11eb-0a83-b1982d8eb146
# ╟─3d96d96e-484d-11eb-0293-b1e8656110fa
# ╟─fa566620-484d-11eb-08ea-cf27826e5d58
# ╟─46523bb0-4850-11eb-2c9b-5d7b5d828246
# ╟─4cc0cf1e-4850-11eb-24e1-6107648b653b
# ╟─450d50e0-4856-11eb-0275-13db9e520e93
# ╠═8b082da0-4391-11eb-1ee3-c70c33b4cf2e
# ╟─1a12ca10-4391-11eb-38f3-475632625048
# ╟─07cf29ee-4b31-11eb-3740-07c6145ffafd
# ╟─4659a0d0-4391-11eb-08bc-5f0fa8380c40
# ╟─65f4fca0-49b3-11eb-2ee4-072e50d481b5
# ╠═79975630-4865-11eb-1c55-290ee27b9677
# ╟─51206ad0-4391-11eb-3d4c-3ba8cbc17a60
# ╠═583db4d2-4391-11eb-2ae9-c73fd33c55a9
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
# ╠═6a63bc30-4865-11eb-0396-6d4ca5f1bb07
# ╟─92d8d6b0-49b3-11eb-2b82-07a4ab8dcaee
@testset "SpatialElem" begin
    @test SI.check_hasmbr(SI.Rect{Float64,2}, SpatialElem{Float64,2,Int,Int})
    @test_throws ArgumentError SI.check_hasmbr(SI.Rect{Float64,2}, Int)
    @test_throws ArgumentError SI.check_hasmbr(SI.Rect{Float64,3}, SpatialElem{Float64,2,Int,Int})
    @test_throws ArgumentError SI.check_hasmbr(SI.Rect{Float64,3}, SpatialElem{Int,2,Int,Int})

    @test SI.check_hasid(Int, SpatialElem{Float64,2,Int,Int})
    @test_throws ArgumentError SI.check_hasid(String, SpatialElem{Float64,2,Int,Int})
    @test_throws ArgumentError SI.check_hasid(Int, SpatialElem{Float64,2,Nothing,Int})
end

@testset "RTree" begin
tree_vars = [SI.RTreeStar, SI.RTreeLinear, SI.RTreeQuadratic]

@testset "Basic Operations" begin
    @testset "RTree{Float,2,Int32,Int}(variant=$tree_var)" for tree_var in tree_vars
        ambr = SI.Rect((0.0, 0.0), (0.0, 0.0))
        bmbr = SI.Rect((0.0, 1.0), (0.0, 1.0))
        cmbr = SI.Rect((0.5, 0.5), (0.5, 0.6))

        tree = RTree{Float64,2}(Int32, Int, variant=tree_var)
        @test tree isa RTree{Float64,2,SpatialElem{Float64,2,Int32,Int}}
        @test SI.variant(tree) === tree_var

        @test eltype(tree) === SpatialElem{Float64,2,Int32,Int}
        @test SI.regiontype(tree) === SI.Rect{Float64,2}
        @test SI.dimtype(tree) === Float64
        @test ndims(tree) === 2
        @test length(tree) == 0
        @test SI.height(tree) == 1
        @test isempty(tree)
        @test isempty(tree, ambr)
        @test isequal(SI.mbr(tree.root), SI.empty(SI.regiontype(tree)))
        @test SI.check(tree)
        @test_throws KeyError delete!(tree, SI.empty(SI.regiontype(tree)), 1)
        @test SI.check(tree)
        @test iterate(tree) === nothing
        @test iterate(contained_in(tree, SI.empty(SI.mbrtype(tree)))) === nothing
        @test iterate(intersects_with(tree, SI.empty(SI.mbrtype(tree)))) === nothing
        @test iterate(contained_in(tree, cmbr)) === nothing
        @test iterate(intersects_with(tree, cmbr)) === nothing
        @test collect(tree) == eltype(tree)[]
        @test collect(contained_in(tree, cmbr)) == eltype(tree)[]
        @test collect(intersects_with(tree, cmbr)) == eltype(tree)[]
        @test typeof(collect(tree)) === Vector{eltype(tree)}

        @test insert!(tree, ambr, 1, 2) === tree
        @test length(tree) == 1
        @test !isempty(tree)
        @test !isempty(tree, ambr)
        @test SI.height(tree) == 1
        @test isequal(SI.mbr(tree.root), ambr)
        @test SI.check(tree)
        @test_throws KeyError delete!(tree, bmbr, 1)
        @test_throws KeyError delete!(tree, ambr, 2)
        @test SI.check(tree)
        @test collect(tree) == [SpatialElem(ambr, Int32(1), 2)]

        @test insert!(tree, bmbr, 2, 2) === tree
        @test length(tree) == 2
        @test SI.height(tree) == 1
        @test isequal(SI.mbr(tree.root), SI.combine(ambr, bmbr))
        @test SI.check(tree)
        @test length(collect(tree)) == 2
        @test_throws KeyError delete!(tree, bmbr, 1)
        @test_throws KeyError delete!(tree, ambr, 2)
        @test_throws KeyError delete!(tree, SI.empty(SI.regiontype(tree)), 1)
        @test delete!(tree, ambr, 1) === tree
        @test length(tree) == 1
        @test SI.check(tree)

        insert!(tree, ambr, 1, 2)
        insert!(tree, cmbr, 2, 3)
        @test length(tree) == 3
        @test_throws SpatialIndexException SI.check(tree) # duplicate ID=2 (b and c)
        @test delete!(tree, cmbr, 2) === tree
        @test length(tree) == 2
        # prepare tree with a(id=1), b(id=2), c(id=3)
        @test insert!(tree, cmbr, 3, 3) === tree
        @test length(tree) == 3
        @test SI.check(tree)

        a = (0.0, 0.0)
        b = (0.55, 0.55)
        c = (0.6, 0.6)
        d = (1.0, 1.0)
        @test length(collect(contained_in(tree, SI.Rect(a, d)))) == 3
        @test length(collect(intersects_with(tree, SI.Rect(a, d)))) == 3
        @test length(collect(contained_in(tree, SI.Rect(a, c)))) == 2
        @test length(collect(intersects_with(tree, SI.Rect(a, c)))) == 2
        @test length(collect(contained_in(tree, SI.Rect(a, b)))) == 1 # a only
        @test length(collect(intersects_with(tree, SI.Rect(a, b)))) == 2 # a and c

        tree2 = similar(tree)
        @test typeof(tree2) === typeof(tree)
        @test tree2 !== tree
        @test isempty(tree2)

        @testset "findfirst()" begin
            @test findfirst(tree, ambr, 2) === nothing
            @test_throws MethodError findfirst(tree, ambr, "1") # wrong key type
            @test_throws MethodError findfirst(tree, SI.empty(SI.Rect{Int,2}()))
            @test findfirst(tree, ambr, 1) == (tree.root, 2)
            @test findfirst(tree, ambr) == (tree.root, 2) # search without id ignores it
            @test findfirst(tree, bmbr, 1) === nothing
            @test findfirst(tree, bmbr, 2) == (tree.root, 1)
            @test findfirst(tree, bmbr) == (tree.root, 1)
            @test findfirst(tree, cmbr, 2) === nothing
            @test findfirst(tree, cmbr, 3) == (tree.root, 3)
            @test findfirst(tree, cmbr) == (tree.root, 3)
            @test findfirst(tree, SI.combine(ambr, cmbr)) === nothing
            @test findfirst(tree, SI.empty(SI.mbrtype(tree))) === nothing
        end
    end

    @testset "RTree{Int,3,String,Nothing}(variant=$tree_var) (no id)" for tree_var in tree_vars
        ambr = SI.Rect((0, 0, 0), (0, 0, 0))
        bmbr = SI.Rect((0, 1, 1), (0, 1, 1))
        tree = RTree{Int,3}(String, variant=tree_var)
        @test tree isa RTree{Int,3,SpatialElem{Int,3,Nothing,String}}
        @test SI.variant(tree) === tree_var
        @test eltype(tree) === SpatialElem{Int,3,Nothing,String}
        @test SI.regiontype(tree) === SI.Rect{Int,3}
        @test SI.dimtype(tree) === Int
        @test ndims(tree) === 3
        @test length(tree) == 0
        @test SI.height(tree) == 1
        @test isempty(tree)
        @test isempty(tree, ambr)
        @test isequal(SI.mbr(tree.root), SI.empty(SI.regiontype(tree)))
        @test SI.check(tree)
        @test_throws KeyError delete!(tree, SI.empty(SI.regiontype(tree)), 1)
        @test_throws KeyError delete!(tree, SI.empty(SI.regiontype(tree)))
        @test SI.check(tree)

        @test_throws Exception insert!(tree, ambr, 1, "2")
        @test insert!(tree, ambr, "2") === tree
        @test length(tree) == 1
        @test !isempty(tree)
        @test !isempty(tree, ambr)
        @test SI.height(tree) == 1
        @test isequal(SI.mbr(tree.root), ambr)
        @test SI.check(tree)
        @test_throws MethodError delete!(tree, bmbr, 1)
        @test_throws KeyError delete!(tree, bmbr)
        @test SI.check(tree)

        @test insert!(tree, bmbr, "3") === tree
        @test length(tree) == 2
        @test SI.height(tree) == 1
        @test isequal(SI.mbr(tree.root), SI.combine(ambr, bmbr))
        @test SI.check(tree)
        @test_throws MethodError delete!(tree, bmbr, 1)
        @test_throws MethodError delete!(tree, SI.empty(SI.Rect{Int,2}))
        @test_throws MethodError delete!(tree, SI.empty(SI.Rect{Float64,3}))
        @test_throws KeyError delete!(tree, SI.empty(SI.regiontype(tree)))
        @test delete!(tree, bmbr) === tree
        @test length(tree) == 1
        @test SI.check(tree)
        @test insert!(tree, bmbr, "4") === tree
        @test length(tree) == 2
        @test SI.check(tree)

        @testset "findfirst()" begin
            @test_throws MethodError findfirst(tree, ambr, 2) # no key
            @test_throws MethodError findfirst(tree, ambr, "1")
            @test_throws MethodError findfirst(tree, SI.empty(SI.Rect{Int,2}()))
            @test_throws MethodError findfirst(tree, SI.empty(SI.Rect{Float64,3}()))
            @test findfirst(tree, ambr) == (tree.root, 1)
            @test findfirst(tree, bmbr) == (tree.root, 2)
            @test findfirst(tree, SI.combine(ambr, bmbr)) === nothing
            @test findfirst(tree, SI.empty(SI.mbrtype(tree))) === nothing
        end
    end
end

@testset "add 1000 vertices" begin
    # generate random point cloud
    Random.seed!(32123)
    mbrs = Vector{SI.Rect{Float64,3}}()
    for i in 1:1000
        w, h, d = 5 .* rand(3)
        x, y, z = 50 .* randn(3)
        rmbr = SI.Rect((x - w, y - h, z - d), (x + w, y + h, z + d))
        push!(mbrs, rmbr)
    end

    @testset "sequential inserts into RTree(variant=$tree_var)" for tree_var in tree_vars
        tree = RTree{SI.dimtype(eltype(mbrs)),ndims(eltype(mbrs))}(Int, String, leaf_capacity=20, branch_capacity=20, variant=tree_var)
        @test tree isa RTree{Float64,3,SpatialElem{Float64,3,Int,String}}
        for (i, rmbr) in enumerate(mbrs)
            insert!(tree, rmbr, i, string(i))
            @test length(tree) == i
            @test tree.nelem_insertions == i
            @test SI.check(tree)
            #@debug "$i: len=$(length(tree)) height=$(SI.height(tree))"
        end
        #@show SI.height(tree)
        #@show tree.nnodes_perlevel

        @testset "findfirst()" begin
            # check that the elements can be found
            for i in 1:length(tree)
                node_ix = findfirst(tree.root, mbrs[i], i)
                @test node_ix !== nothing
                if node_ix !== nothing
                    node, ix = node_ix
                    elem = node[ix]
                    @test SI.id(elem) == i
                    @test elem.val == string(i)
                end
            end
        end

        @testset "iterate" begin
            all_elems = collect(tree)
            @test length(all_elems) == length(tree)
            @test eltype(all_elems) === eltype(tree)

            bound_mbr = SI.Rect((-40.0, -20.0, -20.0), (30.0, 50.0, 40.0))

            in_elems = collect(contained_in(tree, bound_mbr))
            @test eltype(in_elems) === eltype(tree)
            @test length(in_elems) == sum(br -> in(br, bound_mbr), mbrs)

            isect_elems = collect(intersects_with(tree, bound_mbr))
            @test eltype(isect_elems) === eltype(tree)
            @test length(isect_elems) == sum(br -> SI.intersects(br, bound_mbr), mbrs)
        end
    end

    @testset "load!(RTree(variant=$tree_var)) (OMT bulk load)" for tree_var in tree_vars
        tree = RTree{SI.dimtype(eltype(mbrs)),ndims(eltype(mbrs))}(Int, String,
            leaf_capacity=20, branch_capacity=20, variant=tree_var)
        @test tree === SI.load!(tree, enumerate(mbrs), method=:OMT,
            convertel=x -> eltype(tree)(x[2], x[1], string(x[1])))
        @test SI.check(tree)
        @test length(tree) == length(mbrs)
        #@show SI.height(tree)
        #@show tree.nnodes_perlevel
        # cannot bulk-load into non-empty tree
        @test_throws ArgumentError SI.load!(tree, enumerate(mbrs), method=:OMT,
            convertel=x -> eltype(tree)(x[2], x[1], string(x[1])))
        tree2 = similar(tree)
        # can load from tree into another tree
        @test SI.load!(tree2, tree) == tree2
        @test length(tree2) == length(tree)

        @testset "findfirst()" begin
            # check that the elements can be found
            for i in 1:length(tree)
                node_ix = findfirst(tree.root, mbrs[i], i)
                @test node_ix !== nothing
                if node_ix !== nothing
                    node, ix = node_ix
                    elem = node[ix]
                    @test SI.id(elem) == i
                    @test elem.val == string(i)
                end
            end
        end
    end
end

@testset "subtract!(RTree(variant=$tree_var))" for tree_var in tree_vars
    @testset "simple" begin
        tree = RTree{Int,2}(Int, String, leaf_capacity=5, branch_capacity=5, variant=tree_var)
        pts = [(0, 0), (1, 0), (2, 2), (2, 0), (0, 1), (1, 1), (-1, -1)]
        SI.load!(tree, enumerate(pts),
            convertel=x -> eltype(tree)(SI.Rect(x[2], x[2]), x[1], string(x[1])))
        @test length(tree) == length(pts)
        @test SI.check(tree)
        rect = SI.Rect((0, 0), (1, 1))
        @test !isempty(tree, rect)
        @test isempty(tree, SI.Rect((-2, 0), (-1, 1)))
        SI.subtract!(tree, rect)
        @test SI.check(tree)
        @test length(tree) == 3
        @test tree.nelem_deletions == 4
        for i in [3, 4, 7] # check that the correct points stayed in the tree
            @test findfirst(tree.root, SI.Rect(pts[i], pts[i]), i) !== nothing
        end
    end

    @testset "from 1000 points" begin
        # generate random point cloud
        Random.seed!(32123)
        pts = [ntuple(_ -> 5.0 * randn(), 3) for _ in 1:1000]
        tree = RTree{Float64,3}(Int, String, leaf_capacity=10, branch_capacity=10, variant=tree_var)
        SI.load!(tree, enumerate(pts),
            convertel=x -> eltype(tree)(SI.Rect(x[2], x[2]), x[1], string(x[1])))
        @test length(tree) == length(pts)
        #@show tree.nelems tree.nnodes_perlevel

        rect = SI.Rect((-3.0, -4.0, -5.0), (5.0, 8.0, 3.0))
        n_in_rect = sum(pt -> in(SI.Point(pt), rect), pts)
        SI.subtract!(tree, rect)
        @test SI.check(tree)
        @test tree.nelem_deletions == n_in_rect
        #@show tree.nelems tree.nnodes_perlevel tree.nelem_deletions tree.nelem_insertions
    end

    # test how well the R-tree can stand the deletion of all but a single element for various tree sizes
    @testset "subtract!() removing all but single point" begin
        Random.seed!(32123)
        pts = [ntuple(_ -> rand(), 2) for _ in 1:200]
        reftree = RTree{Float64,2}(Int, String, leaf_capacity=5, branch_capacity=5, variant=tree_var)
        loner = eltype(reftree)(SI.Rect((5.0, 5.0), (5.0, 5.0)), 0, "loner")
        @testset "removing $n points" for (n, pt) in enumerate(pts)
            insert!(reftree, eltype(reftree)(SI.Rect(pt, pt), n, string(n)))
            tree = deepcopy(reftree)
            insert!(tree, loner)
            @test length(tree) == n + 1
            SI.subtract!(tree, SI.Rect((0.0, 0.0), (1.0, 1.0)))
            @test length(tree) == 1
            @test SI.check(tree)
            @test SI.id(first(tree)) == 0
        end
    end

    @testset "subtract!() removing n<=N central nodes" begin
        Random.seed!(32123)
        pts = [ntuple(_ -> 2rand() - 1, 3) for _ in 1:200]
        # the test implies that all pts have different distances to the origin
        sort!(pts, by=pt -> maximum(abs, pt))
        reftree = RTree{Float64,3}(Int, String, leaf_capacity=5, branch_capacity=5, variant=tree_var)
        SI.load!(reftree, enumerate(pts),
            convertel=x -> eltype(reftree)(SI.Rect(x[2], x[2]), x[1], string(x[1])))
        corembr = SI.Rect((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))
        tree1 = deepcopy(reftree)
        @test length(tree1) == length(reftree)
        @testset "removing $n points" for (n, pt) in enumerate(pts)
            corembr = SI.combine(corembr, SI.Rect(pt, pt))

            # subtracting centrermbr from tree1 removes just 1 point
            @test length(tree1) == length(reftree) - n + 1
            tree1_to_remove = collect(contained_in(tree1, corembr))
            @test length(tree1_to_remove) == 1
            @test SI.id(first(tree1_to_remove)) == n

            SI.subtract!(tree1, corembr)
            @test length(tree1) == length(reftree) - n
            @test SI.check(tree1)

            # subtracting centrermbr from tree2 removes n points
            tree2 = deepcopy(reftree)
            SI.subtract!(tree2, corembr)
            @test length(tree2) == length(reftree) - n
            @test SI.check(tree2)
        end
    end
end

@testset "show() and print()" begin
    mbrs = Vector{SI.Rect{Float64,2}}(
        [
        SI.Rect((0.0, 0.0), (2.0, 2.0)),
        SI.Rect((-1.0, -1.0), (1.0, 1.0)),
        SI.Rect((-1.0, 0.0), (1.0, 1.0)),
        SI.Rect((0.0, -1.0), (1.0, 1.0)),
        SI.Rect((1.0, 1.0), (2.0, 2.0)),
    ]
    )
    tree = RTree{SI.dimtype(eltype(mbrs)),ndims(eltype(mbrs))}(Int, String, leaf_capacity=4, branch_capacity=4)
    for (i, rmbr) in enumerate(mbrs)
        insert!(tree, rmbr, i, string(i))
    end

    eltyp = eltype(tree)
    rectyp = eltype(mbrs)
    test_show_string = "$(typeof(tree))(variant=RTreeStar, tight_mbrs=true, nearmin_overlap=1, fill_factor=0.7, split_factor=0.4, reinsert_factor=0.3, leaf_capacity=4, branch_capacity=4)\n5 element(s) in 2 level(s) (1, 2 node(s) per level):\n level=2 nchildren=2 mbr=((-1.0, -1.0), (2.0, 2.0))"
    test_print_string = "$(typeof(tree))(variant=RTreeStar, tight_mbrs=true, nearmin_overlap=1, fill_factor=0.7, split_factor=0.4, reinsert_factor=0.3, leaf_capacity=4, branch_capacity=4)\n5 element(s) in 2 level(s) (1, 2 node(s) per level):\n level=2 nchildren=2 mbr=((-1.0, -1.0), (2.0, 2.0)):\n  level=1 nchildren=3 mbr=((-1.0, -1.0), (1.0, 1.0)):\n   $eltyp($rectyp((-1.0, -1.0), (1.0, 1.0)), 2, \"2\")\n   $eltyp($rectyp((-1.0, 0.0), (1.0, 1.0)), 3, \"3\")\n   $eltyp($rectyp((0.0, -1.0), (1.0, 1.0)), 4, \"4\")\n  level=1 nchildren=2 mbr=((0.0, 0.0), (2.0, 2.0)):\n   $eltyp($rectyp((0.0, 0.0), (2.0, 2.0)), 1, \"1\")\n   $eltyp($rectyp((1.0, 1.0), (2.0, 2.0)), 5, \"5\")"

    io = IOBuffer()
    show(io, tree)
    @test String(take!(io)) == test_show_string
    print(io, tree)
    @test String(take!(io)) == test_print_string
    show(io, tree; recurse=false)
    @test String(take!(io)) == test_show_string
    show(io, tree; recurse=true)
    @test String(take!(io)) == test_print_string
end
end
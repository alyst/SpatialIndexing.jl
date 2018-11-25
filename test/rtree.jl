@testset "RTree" begin
@testset "SpatialElem" begin
    @test SI.check_hasmbr(SI.Rect{Float64, 2}, SpatialElem{Float64, 2, Int, Int})
    @test_throws ArgumentError SI.check_hasmbr(SI.Rect{Float64, 2}, Int)
    @test_throws ArgumentError SI.check_hasmbr(SI.Rect{Float64, 3}, SpatialElem{Float64, 2, Int, Int})
    @test_throws ArgumentError SI.check_hasmbr(SI.Rect{Float64, 3}, SpatialElem{Int, 2, Int, Int})

    @test SI.check_hasid(Int, SpatialElem{Float64, 2, Int, Int})
    @test_throws ArgumentError SI.check_hasid(String, SpatialElem{Float64, 2, Int, Int})
    @test_throws ArgumentError SI.check_hasid(Int, SpatialElem{Float64, 2, Nothing, Int})
end

@testset "Basic Operations" begin
    tree_vars = [SI.RTreeStar, SI.RTreeLinear, SI.RTreeQuadratic]
    @testset "RTree{Float,2,Int32,Int}(variant=$tree_var)" for tree_var in tree_vars
        ambr = SI.Rect((0.0, 0.0), (0.0, 0.0))
        bmbr = SI.Rect((0.0, 1.0), (0.0, 1.0))

        tree = RTree{Float64, 2}(Int32, Int, variant=tree_var)
        @test tree isa RTree{Float64, 2, SpatialElem{Float64, 2, Int32, Int}}
        @test SI.variant(tree) === tree_var

        @test eltype(tree) === SpatialElem{Float64, 2, Int32, Int}
        @test SI.regiontype(tree) === SI.Rect{Float64, 2}
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

        @test insert!(tree, bmbr, 2, 2) === tree
        @test length(tree) == 2
        @test SI.height(tree) == 1
        @test isequal(SI.mbr(tree.root), SI.combine(ambr, bmbr))
        @test SI.check(tree)
        @test_throws KeyError delete!(tree, bmbr, 1)
        @test_throws KeyError delete!(tree, ambr, 2)
        @test_throws KeyError delete!(tree, SI.empty(SI.regiontype(tree)), 1)
        @test delete!(tree, ambr, 1) === tree
        @test length(tree) == 1
        @test SI.check(tree)
    end

    @testset "RTree{Int,3,String,Nothing}(variant=$tree_var) (no id)" for tree_var in tree_vars
        ambr = SI.Rect((0, 0, 0), (0, 0, 0))
        bmbr = SI.Rect((0, 1, 1), (0, 1, 1))
        tree = RTree{Int, 3}(String, variant=tree_var)
        @test tree isa RTree{Int, 3, SpatialElem{Int, 3, Nothing, String}}
        @test SI.variant(tree) === tree_var
        @test eltype(tree) === SpatialElem{Int, 3, Nothing, String}
        @test SI.regiontype(tree) === SI.Rect{Int, 3}
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

        @test_throws MethodError insert!(tree, ambr, 1, "2")
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
        @test_throws MethodError delete!(tree, SI.empty(SI.Rect{Int, 2}))
        @test_throws MethodError delete!(tree, SI.empty(SI.Rect{Float64, 3}))
        @test_throws KeyError delete!(tree, SI.empty(SI.regiontype(tree)))
        @test delete!(tree, bmbr) === tree
        @test length(tree) == 1
        @test SI.check(tree)
    end
end

@testset "1000 vertices" begin
    # generate random point cloud
    Random.seed!(32123)
    mbrs = Vector{SI.Rect{Float64, 3}}()
    for i in 1:1000
        w, h, d = 5 .* rand(3)
        x, y, z = 50 .* randn(3)
        rmbr = SI.Rect((x-w, y-h, z-d), (x+w, y+h, z+d))
        push!(mbrs, rmbr)
    end

    @testset "sequential inserts" begin
        tree = RTree{SI.dimtype(eltype(mbrs)), ndims(eltype(mbrs))}(Int, String, leaf_capacity = 20, branch_capacity = 20)
        @test tree isa RTree{Float64, 3, SpatialElem{Float64, 3, Int, String}}
        for (i, rmbr) in enumerate(mbrs)
            insert!(tree, rmbr, i, string(i))
            @test length(tree) == i
            @test tree.nelem_insertions == i
            @test SI.check(tree)
            #@debug "$i: len=$(length(tree)) height=$(SI.height(tree))"
        end
        @show SI.height(tree)
        @show tree.nnodes_perlevel

        @testset "findleaf" begin
            # check that the elements can be found
            for i in 1:length(tree)
                node_ix = SI.findleaf(tree.root, mbrs[i], i)
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

    @testset "bulk load" begin
        tree = RTree{SI.dimtype(eltype(mbrs)), ndims(eltype(mbrs))}(Int, String,
                                leaf_capacity = 20, branch_capacity = 20)
        @test tree === SI.load!(tree, enumerate(mbrs), method=:OMT,
                                getid = x -> x[1], getmbr = x -> x[2], getval = x -> string(x[1]))
        @test SI.check(tree)
        @test length(tree) ==length(mbrs)
        @show SI.height(tree)
        @show tree.nnodes_perlevel
        # cannot bulk-load into non-empty tree
        @test_throws ArgumentError SI.load!(tree, enumerate(mbrs), method=:OMT,
                                            getid = x -> x[1], getmbr = x -> x[2], getval = x -> string(x[1]))

        @testset "findleaf" begin
            # check that the elements can be found
            for i in 1:length(tree)
                node_ix = SI.findleaf(tree.root, mbrs[i], i)
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

@testset "subtract!()" begin
    @testset "simple" begin
        tree = RTree{Int, 2}(Int, String, leaf_capacity = 5, branch_capacity = 5)
        pts = [(0, 0), (1, 0), (2, 2), (2, 0), (0, 1), (1, 1), (-1, -1)]
        SI.load!(tree, enumerate(pts),
                 getid = x -> x[1],
                 getmbr = x -> SI.Rect(x[2], x[2]),
                 getval = x -> string(x[1]))
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
            @test SI.findleaf(tree.root, SI.Rect(pts[i], pts[i]), i) !== nothing
        end
    end

    @testset "from 1000 points" begin
        # generate random point cloud
        Random.seed!(32123)
        pts = [ntuple(_ -> 5. * randn(), 3) for _ in 1:1000]
        tree = RTree{Float64, 3}(Int, String, leaf_capacity = 10, branch_capacity = 10)
        SI.load!(tree, enumerate(pts),
                 getid = x -> x[1], getmbr = x -> SI.Rect(x[2], x[2]), getval = x -> string(x[1]))
        @test length(tree) == length(pts)
        @show tree.nelems tree.nnodes_perlevel

        rect = SI.Rect((-3.0, -4.0, -5.0), (5.0, 8.0, 3.0))
        n_in_rect = sum(pt -> in(SI.Point(pt), rect), pts)
        SI.subtract!(tree, rect)
        @test SI.check(tree)
        @test tree.nelem_deletions == n_in_rect
        @show tree.nelems tree.nnodes_perlevel tree.nelem_deletions tree.nelem_insertions
    end
end

end

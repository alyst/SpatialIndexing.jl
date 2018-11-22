@testset "Regions" begin

@testset "Point" begin
    @test SI.Point((1, 2)) isa SI.Point{Int, 2}
    @test SI.Point((1.0, 2.0, 3.0)) isa SI.Point{Float64, 3}
    @test SI.Point{Float64, 3}((1, 2, 3)) isa SI.Point{Float64, 3}
    @test isequal(SI.empty(SI.Point{Float64, 2}).coord, (NaN, NaN))
    @test isequal(SI.empty(SI.Point{Int, 1}).coord, (0,))
    @test SI.isvalid(SI.Point((1, 2)))
    @test SI.isvalid(SI.Point((1.0, 2.0, 3.0)))
    @test !SI.isvalid(SI.Point((NaN, 1.0)))

    @test_throws MethodError SI.Point((1.0, 2, 3.0))
    @test_throws InexactError SI.Point{Int,2}((1.5, 3.0))
    @test_throws MethodError SI.Point{Int,3}((1, 3))

    @test SI.area(SI.Point((1,2))) === 0
    @test SI.area(SI.Point((1.0,2.0))) === 0.0
    @test SI.perimeter(SI.Point((1,2))) === 0
    @test SI.perimeter(SI.Point((1.0,2.0))) === 0.0
    @test SI.sqrdistance(SI.Point((1,)), SI.Point((3,))) === 4
    @test SI.sqrdistance(SI.Point((1, 0)), SI.Point((3, -3))) === 13
end

@testset "Rect" begin
    @test SI.Rect((1, 2), (2, 3)) isa SI.Rect{Int, 2}
    @test SI.Rect((1.0, 2.0, 3.0), (2, 3, 4)) isa SI.Rect{Float64, 3}
    @test SI.Rect{Float64, 3}((1, 2, 3), (2, 3, 4)) isa SI.Rect{Float64, 3}
    @test isequal(SI.empty(SI.Rect{Float64, 2}), SI.Rect((NaN, NaN), (NaN, NaN)))
    @test isequal(SI.empty(SI.Rect{Int, 1}), SI.Rect((0,), (0,)))
    @test SI.isvalid(SI.Rect((1, 2), (2, 3)))
    @test SI.isvalid(SI.Rect((1.0, 2.0), (2.0, 3.0)))
    @test SI.isvalid(SI.Rect((1.0, 2.0), (Inf, Inf)))
    @test SI.isvalid(SI.Rect((-Inf, -Inf), (2.0, 3.0)))
    @test !SI.isvalid(SI.Rect((NaN, 2.0), (2.0, 3.0)))
    @test !SI.isvalid(SI.Rect((1.0, 2.0), (2.0, NaN)))
    @test !SI.isvalid(SI.Rect((2.5, 2.0), (2.0, 3.0)))
    @test !SI.isvalid(SI.Rect((1.0, 2.0), (2.0, -3.0)))
    @test !SI.isvalid(SI.Rect((1.0, Inf), (2.0, 3.0)))

    @test_throws MethodError SI.Rect((1.0, 2, 3.0), (2, 3.0, 4))
    @test_throws MethodError SI.Rect((1.0, 3.0), (2, 3.0, 4))
    @test_throws InexactError SI.Rect{Int,2}((1.0, 3.0), (1.5, 2.0))

    @test SI.area(SI.Rect((1.0,), (2.0,))) === 1.0
    @test SI.area(SI.Rect((1,), (0,))) === -1
    @test SI.area(SI.Rect((1.0, 2.0), (2.5, 4.0))) === 3.0

    @test SI.perimeter(SI.Rect((1.0,), (2.0,))) === 1.0
    @test SI.perimeter(SI.Rect((1,), (0,))) === -1
    @test SI.perimeter(SI.Rect((1.0, 2.0), (2.5, 4.0))) === 7.0

    @test SI.center(SI.Rect((1.0,), (2.0,))) == SI.Point((1.5,))
    @test SI.center(SI.Rect((1.0, 2.0), (2.5, 4.0))) == SI.Point((1.75, 3.0))

    ai = SI.Rect((0, 0), (1, 1))
    bi = SI.Rect((1, 1), (2, 2))
    ci = SI.Rect((2, 2), (3, 3))
    @test SI.intersect(ai, bi) == SI.Rect((1,1), (1,1))
    @test SI.intersect(ai, ci) === SI.empty(SI.Rect{Int, 2})
    @test SI.combine(ai, bi) == SI.Rect((0,0), (2,2))
    @test SI.center(SI.Rect((2, 2), (3, 4))) == SI.Point((2.5, 3.0))
    @test SI.center(SI.Rect((1, 2), (3, 4))) == SI.Point((2.0, 3.0))

    a = SI.Rect{Float64, 2}((0, 0), (1, 1))
    b = SI.Rect{Float64, 2}((1, 1), (2, 2))
    c = SI.Rect{Float64, 2}((2, 2), (3, 3))
    @test SI.intersect(a, b) == SI.Rect((1.0, 1.0), (1.0, 1.0))
    @test SI.intersect(a, c) === SI.empty(typeof(a))
    @test SI.combine(a, b) == SI.Rect((0.0,0.0), (2.0,2.0))
    @test SI.combined_area(a, b) == 4
    @test SI.combined_area(a, c) == 9
    @test SI.overlap_area(a, b) == 0
    @test SI.overlap_area(a, c) == 0
    @test SI.overlap_area(SI.combine(a, b), SI.combine(b, c)) == 1

    @test SI.touches(a, SI.combine(a, b))
    @test SI.touches(SI.combine(a, b), a)
    @test SI.touches(b, SI.combine(a, b))
    @test !SI.touches(SI.combine(a, b), c)
    @test SI.touches(SI.combine(a, c), a)
    @test SI.touches(SI.combine(a, c), c)
    @test !SI.touches(SI.combine(a, c), b)
    @test SI.touches(SI.combine(a, b), SI.Rect((1.0, 0.0), (2.0, 1.0)))
    @test SI.touches(SI.combine(a, b), SI.Rect((-1.0, 0.0), (0.0, 1.0))) #outside
    @test SI.touches(SI.combine(a, b), SI.Rect((0.0, -1.0), (1.0, 0.0))) # outside
    @test SI.touches(SI.combine(a, b), SI.Rect((2.0, 0.0), (3.0, 1.0))) # outside
    @test SI.touches(SI.combine(a, b), SI.Rect((0.0, 2.0), (1.0, 3.0))) # outside
end

end

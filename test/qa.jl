using StatsFuns
using Aqua: Aqua
using ExplicitImports: ExplicitImports
using JET: JET
using Test

@testset "QA" begin
    @testset "Aqua" begin
        Aqua.test_all(StatsFuns)
    end

    @testset "ExplicitImports" begin
        # No implicit imports (`using XY`)
        @test ExplicitImports.check_no_implicit_imports(StatsFuns) === nothing

        # All explicit imports (`using XY: Z`) are loaded via their owners
        @test ExplicitImports.check_all_explicit_imports_via_owners(
            StatsFuns;
            ignore = (
                # Ref https://github.com/JuliaTesting/ExplicitImports.jl/issues/92
                :digamma,
            ),
        ) === nothing

        # Limit explicit imports (`using XY: Z`) of non-public names to a minimum
        @test ExplicitImports.check_all_explicit_imports_are_public(
            StatsFuns;
            ignore = (
                # Ref https://github.com/JuliaTesting/ExplicitImports.jl/issues/92
                :digamma,
            ),
        ) === nothing

        # No explicit imports (`using XY: Z`) that are not used
        @test ExplicitImports.check_no_stale_explicit_imports(StatsFuns; ignore = (:digamma,)) === nothing

        # Nothing is accessed via modules other than its owner
        @test ExplicitImports.check_all_qualified_accesses_via_owners(StatsFuns) === nothing

        # Limit accesses of non-public names to a minimum
        @test ExplicitImports.check_all_qualified_accesses_are_public(
            StatsFuns;
            ignore = (
                (VERSION < v"1.11" ? (:Fix2,) : ())...,
                :promote_typeof,
            ),
        ) === nothing

        # No self-qualified accesses
        @test ExplicitImports.check_no_self_qualified_accesses(StatsFuns) === nothing
    end

    @testset "JET" begin
        # Check that there are no undefined global references and undefined field accesses
        JET.test_package(StatsFuns; target_defined_modules = true, mode = :typo)

        # Analyze methods based on their declared signature
        JET.report_package(StatsFuns; target_defined_modules = true)
    end
end

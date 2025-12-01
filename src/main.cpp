import std;
import types;
import math;
import muller;
import krylov;


#define EXAMPLE 1

auto example() -> void {
    #if EXAMPLE == 1
    std::vector<std::complex<f64>> dataA = {
        { 5, 0 }, { 3, 0 }, { 1, 0 }, { 3, 0 },
        { 2, 0 }, { 3, 0 }, {-2, 0 }, { 0, 0 },
        {-1, 0 }, { 0, 0 }, { 4, 0 }, { 0, 0 },
        {-1, 0 }, { 0, 0 }, { 1, 0 }, { 3, 0 },
    };
    auto A = Matrix::From(std::move(dataA));
    #elif EXAMPLE == 2
    std::vector<std::complex<f64>> dataA = {
        { 1, 0 }, { 2, 0 }, { 3, 0 }, { 4, 0 },
        { 3, 0 }, { 4, 0 }, { 5, 0 }, { 3, 0 },
        { 1, 0 }, { 1, 0 }, { 4, 0 }, { 1, 0 },
        { 4, 0 }, { 1, 0 }, { 2, 0 }, { 2, 0 },
    };
    auto A = Matrix::From(std::move(dataA));
    #elif EXAMPLE == 3
    auto A = Matrix::Hilbert(6);
    #endif

    std::println("A");
    A.print();

    f64 cond_inf = A.map().rowwise().lpNorm<1>().maxCoeff() * A.map().inverse().rowwise().lpNorm<1>().maxCoeff();
    std::println("cond_inf = {}", cond_inf);

    auto phi = findPoly(A, Vector::Basis(A.size(), 0));
    std::println("phi");
    phi.print();

    auto roots = findRoots(phi);
    std::println("roots");
    for (u64 i = 0; i < roots.size(); ++i) {
        printComplex(roots[i]);
        std::println();
    }

    std::println("missing roots: {}", A.size() - phi.degree());

    auto sum = std::reduce(roots.begin(), roots.end());
    std::println("err");
    std::println("{:10.15e}", std::abs(sum - A.trace()));
}

auto experiment1() -> void {
    u64 n = 10, exp = 100;
    for (u64 k = 0; k < exp; ++k)  {
        auto A = Matrix::Random(n);
        auto NA = std::abs(std::sqrt(A.map().cwiseSquare().sum()));
        auto A2 = std::abs(std::sqrt((A.map().transpose() * A.map()).eigenvalues()[n-1]));
        if (!(NA / std::sqrt(n) <= A2 && A2 <= NA)) {
            std::println("doesn't work");
        }
    }
    std::println("if you saw nothing then you are good");
}

auto experiment2() -> void {
    u64 n = 10, exp = 100;
    for (u64 k = 0; k < exp; ++k)  {
        auto A = Matrix::Random(n);
        auto A22 = std::abs((A.map().transpose() * A.map()).eigenvalues()[n-1]);
        auto A1 = A.map().colwise().lpNorm<1>().maxCoeff();
        auto AINF = A.map().rowwise().lpNorm<1>().maxCoeff();
        if (!(A22 <= A1 * AINF)) {
            std::println("doesn't work");
        }
    }
    std::println("if you saw nothing then you are good");
}

auto main() -> i32 {
    // uncomment what you need

    example();
    // experiment1();
    // experiment2();

    return 0;
}
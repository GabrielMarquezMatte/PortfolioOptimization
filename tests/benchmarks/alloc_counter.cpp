#include "alloc_counter.hpp"

#include <atomic>
#include <cstdio>
#include <cstdlib>
#include <new>
#include <vector>

namespace po::bench
{

namespace
{

std::atomic<std::uint64_t> allocations{0};

void* volatile sink = nullptr;

void escape(const void* pointer)
{
    sink = const_cast<void*>(pointer);
}

int expect_delta(const char* what, std::uint64_t before, std::uint64_t expected)
{
    std::uint64_t delta = allocation_count() - before;
    if (delta == expected)
    {
        return 0;
    }
    std::printf("allocation selftest: %s counted %llu allocations, expected %llu\n", what,
                static_cast<unsigned long long>(delta), static_cast<unsigned long long>(expected));
    return 1;
}

int count_vector_allocations()
{
    std::uint64_t before = allocation_count();
    std::vector<double> values(100);
    escape(values.data());
    int failures = expect_delta("one vector", before, 1);
    before = allocation_count();
    std::vector<double> first(10);
    std::vector<double> second(10);
    std::vector<double> third(10);
    escape(first.data());
    escape(second.data());
    escape(third.data());
    failures += expect_delta("three vectors", before, 3);
    return failures;
}

int count_array_and_nothrow_allocations()
{
    std::uint64_t before = allocation_count();
    double* array = new double[16];
    escape(array);
    delete[] array;
    int failures = expect_delta("new[]", before, 1);
    before = allocation_count();
    int* single = new (std::nothrow) int(7);
    escape(single);
    delete single;
    failures += expect_delta("nothrow new", before, 1);
    return failures;
}

int count_no_allocations()
{
    std::uint64_t before = allocation_count();
    double sum = 0.0;
    double stack[64] = {};
    for (double& value : stack)
    {
        value += 1.0;
        sum += value;
    }
    escape(&sum);
    return expect_delta("stack work", before, 0);
}

} // namespace

std::uint64_t allocation_count()
{
    return allocations.load(std::memory_order_relaxed);
}

int allocation_selftest()
{
    int failures = count_vector_allocations() + count_array_and_nothrow_allocations() + count_no_allocations();
    if (failures == 0)
    {
        std::printf("allocation selftest ok\n");
    }
    return failures;
}

} // namespace po::bench

// GCC cannot see that the replaced operator new above is the one that hands out the pointers freed here.
#if defined(__GNUC__) && !defined(__clang__)
#pragma GCC diagnostic ignored "-Wmismatched-new-delete"
#endif

// The global allocation functions below count every plain allocation of the process. Aligned (over-aligned
// type) allocations keep the library's own functions and are not counted; nothing measured here uses them.
void* operator new(std::size_t size)
{
    po::bench::allocations.fetch_add(1, std::memory_order_relaxed);
    void* pointer = std::malloc(size == 0 ? 1 : size);
    if (pointer == nullptr)
    {
        throw std::bad_alloc();
    }
    return pointer;
}

void* operator new(std::size_t size, const std::nothrow_t&) noexcept
{
    po::bench::allocations.fetch_add(1, std::memory_order_relaxed);
    return std::malloc(size == 0 ? 1 : size);
}

void* operator new[](std::size_t size)
{
    return operator new(size);
}

void* operator new[](std::size_t size, const std::nothrow_t& tag) noexcept
{
    return operator new(size, tag);
}

void operator delete(void* pointer) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, const std::nothrow_t&) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, const std::nothrow_t&) noexcept
{
    std::free(pointer);
}

void operator delete(void* pointer, std::size_t) noexcept
{
    std::free(pointer);
}

void operator delete[](void* pointer, std::size_t) noexcept
{
    std::free(pointer);
}

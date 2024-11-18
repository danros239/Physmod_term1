#include <cmath>
#include <iostream>
#pragma once

// template<typename T>
// class State
// {
//     public:
//     virtual State operator+ (State other);
//     virtual State operator* (T coeff);
// };

template <typename T>
class SingleState // : public State<T>
{
    public:
    T x, v, t;
    SingleState(T x0, T v0): x(x0), v(v0), t(0) { };
    SingleState(T x0, T v0, T t0): x(x0), v(v0), t(t0) { };
    SingleState(): x(0), v(0), t(0) { };

    SingleState operator+ (SingleState other)
    {
        SingleState res = *this;
        res.x += other.x;
        res.v += other.v;
        return res;
    }

    SingleState operator* (T coeff)
    {
        SingleState res;
        res.x = x * coeff;
        res.v = v * coeff;
        return res;
    }

    void print()
    {
        std::cout << x << " " << v << std::endl;
    }
};

template <typename T, template <typename F_T> class Function>
class Abstract_Solver
{
    public:
    void rk_solve(SingleState<T>& state, Function<T> foo, T dt)
    {
        SingleState<T> k1(0,0,state.t), k2(0,0,state.t), k3(0,0,state.t), k4(0,0,state.t);
        k1 = foo(state);
        k2 = foo(state + k1*dt*static_cast<T>(0.5));
        k3 = foo(state + k2*dt*static_cast<T>(0.5));
        k4 = foo(state + k3*dt);
        state = state + (k1*(T)(1/(T)6) + k2*(T)(1/(T)3) + k3*(T)(1/(T)3) + k4*(T)(1/(T)6))*dt;
        state.t += dt;
    }
    Abstract_Solver() { };
};

// template <typename T>
// class Function
// {
//     virtual State<T> operator()(State<T> state);
// };

template<typename T>
class MPendulum // : public Function<T>
{
    public:
    T w;
    SingleState<T> operator()(SingleState<T> state)
    {
        SingleState<T> newstate;
        newstate.v = -state.x*w*w;
        newstate.x = state.v;
        newstate.t = state.t;
        return newstate;
    }

    MPendulum(T omega): w(omega) { };
};
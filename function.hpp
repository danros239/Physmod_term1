#include <cmath>
#include <iostream>
#include <vector>
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

    SingleState ts(T dt) // returns state shifted in time by dt
    {
        SingleState other = *this;
        other.t += dt;
        return other;
    }

    void print()
    {
        std::cout << x << " " << v << std::endl;
    }
    T energy()
    {
        return x*x + v*v;
    }
};

// template<typename T>
// class NState
// {
//     public:
//     const size_t n; // number of coordinates

//     std::vector<T> x, v;
    
//     NState(size_t dim, std::vector x0, std::vector v_0, T t0): n(dim), x(x0), v(v0), t(t0)
//     {

//     }
// };

enum algos
{
    EULER,
    HEUN,
    RK4
};

template <typename T, template <typename F_T> class Function>
class Abstract_Solver
{
    public:
    algos type;
    void euler_solve(SingleState<T>& state, Function<T> foo, T dt)
    {
        SingleState<T> delta = foo(state)*dt;
        state = state + delta;
        state.t += dt;
    }

    void heun_solve(SingleState<T>& state, Function<T> foo, T dt)
    {
        SingleState<T> k1 = foo(state)*dt*(T)0.5;
        SingleState<T> k2 = foo(state.ts(dt*(T)0.5) + k1);
        state = state + k2*(T)dt;
        state.t += dt;
    }

    void rk_solve(SingleState<T>& state, Function<T> foo, T dt)
    {
        SingleState<T> k1(0,0,state.t), k2(0,0,state.t + 0.5*dt), k3(0,0,state.t + 0.5*dt), k4(0,0,state.t*dt);
        k1 = foo(state);
        k2 = foo(state.ts(dt*(T)0.5) + k1*dt*static_cast<T>(0.5));
        k3 = foo(state.ts(dt*(T)0.5) + k2*dt*static_cast<T>(0.5));
        k4 = foo(state.ts(dt) + k3*dt);
        state = state + (k1*(T)(1/(T)6) + k2*(T)(1/(T)3) + k3*(T)(1/(T)3) + k4*(T)(1/(T)6))*dt;
        state.t += dt;
    }
    T total_energy(SingleState<T> state)
    {

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

template <typename T>
class MP_Friction
{
    public:
    T w;
    T alpha;
    SingleState<T> operator()(SingleState<T> state)
    {
        SingleState<T> newstate;
        newstate.v = -state.x*w*w - state.v * alpha;
        newstate.x = state.v;
        newstate.t = state.t;
        return newstate;
    }
    MP_Friction(T omega, T alpha): w(omega), alpha(alpha) { };
};

template <typename T>
class MP_Periodic_Force
{
    public:
    T w0;
    T alpha;
    T f, w;
    SingleState<T> operator()(SingleState<T> state)
    {
        SingleState<T> newstate;
        newstate.v = -state.x*w0*w0 - state.v * alpha + f*cos(w*state.t);
        newstate.x = state.v;
        newstate.t = state.t;
        return newstate;
    }
    MP_Periodic_Force(T omega_0, T omega, T f, T alpha): w(omega), alpha(alpha), f(f), w0(omega_0) { };
};

template <typename T>
class Double_Pendulum
{
    T l1, l2;
    T m1, m2;
    T g;

};
#include <cmath>
#include <iostream>
#include <vector>
#pragma once

template<typename T>
class State
{
    // public:
    // virtual State operator+ (State other);
    // virtual State operator* (T coeff);
};

template <typename T>
class DoubleState : public State<T>
{
    public:
    T x, v, t;
    DoubleState(T x0, T v0): x(x0), v(v0), t(0) { };
    DoubleState(T x0, T v0, T t0): x(x0), v(v0), t(t0) { };
    DoubleState(): x(0), v(0), t(0) { };

    DoubleState operator+ (DoubleState other)
    {
        DoubleState res = *this;
        res.x += other.x;
        res.v += other.v;
        return res;
    }

    DoubleState operator* (T coeff)
    {
        DoubleState res;
        res.x = x * coeff;
        res.v = v * coeff;
        return res;
    }

    DoubleState ts(T dt) // returns state shifted in time by dt
    {
        DoubleState other = *this;
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

template<typename T>
class NState // state of a system with N coordinates
{
    public:
    const size_t n; // number of coordinates

    std::vector<T> x; // a vector of position coordinates
    T t;
    
    NState(size_t dim, std::vector<T> x_0, T t_0): n(dim), x(x_0), t(t_0){ };

    NState operator+(NState other)
    {
        NState res = *this;
        res.x += other.x;
        return res;
    }

    NState operator*(T coeff)
    {
        NState res = *this;
        res.x *= coeff;
        return res;
    }
};

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
    void euler_solve(DoubleState<T>& state, Function<T> foo, T dt)
    {
        DoubleState<T> delta = foo(state)*dt;
        state = state + delta;
        state.t += dt;
    }

    void heun_solve(DoubleState<T>& state, Function<T> foo, T dt)
    {
        DoubleState<T> k1 = foo(state)*dt*(T)0.5;
        DoubleState<T> k2 = foo(state.ts(dt*(T)0.5) + k1);
        state = state + k2*(T)dt;
        state.t += dt;
    }

    void rk_solve(DoubleState<T>& state, Function<T> foo, T dt)
    {
        DoubleState<T> k1(0,0,state.t), k2(0,0,state.t + 0.5*dt), k3(0,0,state.t + 0.5*dt), k4(0,0,state.t*dt);
        k1 = foo(state);
        k2 = foo(state.ts(dt*(T)0.5) + k1*dt*static_cast<T>(0.5));
        k3 = foo(state.ts(dt*(T)0.5) + k2*dt*static_cast<T>(0.5));
        k4 = foo(state.ts(dt) + k3*dt);
        state = state + (k1*(T)(1/(T)6) + k2*(T)(1/(T)3) + k3*(T)(1/(T)3) + k4*(T)(1/(T)6))*dt;
        state.t += dt;
    }
    T total_energy(DoubleState<T> state)
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
    const size_t dim = 2;
    // DoubleState<T> operator()(DoubleState<T> state)
    // {
    //     DoubleState<T> newstate;
    //     newstate.v = -state.x*w*w;
    //     newstate.x = state.v;
    //     newstate.t = state.t;
    //     return newstate;
    // }

    NState<T> operator()(NState<T> state) // calculates time derivative of a state
    {
        NState<T> newstate;
        newstate.x[1] = -state.x[0]*w*w;
        newstate.x[0] = state.x[1];
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
    DoubleState<T> operator()(DoubleState<T> state)
    {
        DoubleState<T> newstate;
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
    DoubleState<T> operator()(DoubleState<T> state)
    {
        DoubleState<T> newstate;
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

template <typename T, template <typename F_T> class Function>
class Algebraic_Solver
{
    T newton_solve()
    {

    }
};
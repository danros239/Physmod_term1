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

template <typename T> class NState
{
    public:

    size_t n;
    std::vector<T> x; // a vector of position coordinates
    T t;
    
    NState(size_t dim, std::vector<T> x_0, T t_0): n(dim), x(x_0), t(t_0)
    {
        x.resize(n);
    };

    NState(size_t dim): n(dim) 
    { 
        x.resize(n);
    };

    NState operator+(NState other)
    {
        NState res = *this;
        for(int i=0; i<n; i++)
            res.x[i] += other.x[i];
        return res;
    }
    NState operator*(T coeff)
    {
        NState res = *this;
        for(int i=0; i<n; i++)
            res.x[i] *= coeff;
        return res;
    }

    NState ts(T timeshift)
    {
        NState other = *this;
        other.t += timeshift;
        return other;
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
    void euler_solve(NState<T>& state, Function<T> foo, T dt)
    {
        NState<T> delta = foo(state)*dt;
        state = state + delta;
        state.t += dt;
    }

    void heun_solve(NState<T>& state, Function<T> foo, T dt)
    {
        NState<T> k1 = foo(state)*dt*(T)0.5;
        NState<T> k2 = foo(state.ts(dt*(T)0.5) + k1);
        state = state + k2*(T)dt;
        state.t += dt;
    }

    void rk_solve(NState<T>& state, Function<T> foo, T dt)
    {
        std::vector<T> tmp{0.0, 0.0};
        NState<T> k1(2,tmp,state.t), k2(2,tmp,state.t + 0.5*dt), k3(2,tmp,state.t + 0.5*dt), k4(2,tmp,state.t*dt);
        k1 = foo(state);
        k2 = foo(state.ts(dt*(T)0.5) + k1*dt*static_cast<T>(0.5));
        k3 = foo(state.ts(dt*(T)0.5) + k2*dt*static_cast<T>(0.5));
        k4 = foo(state.ts(dt) + k3*dt);
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
    const size_t dim = 2;
    // NState<T> operator()(DoubleState<T> state)
    // {
    //     DoubleState<T> newstate;
    //     newstate.v = -state.x*w*w;
    //     newstate.x = state.v;
    //     newstate.t = state.t;
    //     return newstate;
    // }

    NState<T> operator()(NState<T> state) // calculates time derivative of a state
    {
        NState<T> newstate(2);
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
    NState<T> operator()(NState<T> state)
    {
        NState<T> newstate(2);
        newstate.x[1] = -state.x[0]*w*w - state.x[1] * alpha;
        newstate.x[0] = state.x[1];
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
    NState<T> operator()(NState<T> state)
    {
        NState<T> newstate(2);
        newstate.x[1] = -state.x[0]*w0*w0 - state.x[1] * alpha + f*cos(w*state.t);
        newstate.x[0] = state.x[1];
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
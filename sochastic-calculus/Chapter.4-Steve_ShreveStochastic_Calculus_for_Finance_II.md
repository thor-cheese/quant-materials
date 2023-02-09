

\section{Stochastic Calculus}

\subsection{Introduction}

This chapter defines Itô integrals and develops their properties. These are used to model the value of a portfolio that results from trading assets in continuous time. The calculus used to manipulate these integrals is based on the Itô-Doeblin formula of Section $4.4$ and differs from ordinary calculus. This difference can be traced to the fact that Brownian motion has a nonzero quadratic variation and is the source of the volatility term in the Black-Scholes-Merton partial differential equation. The Black-Scholes-Merton equation is presented in Section 4.5. This is in the spirit of Sections $1.1$ and $1.2$ of Volume I in which we priced options by determining the portfolio that would hedge a short position. In particular, there is no discussion of riskneutral pricing in this chapter. That topic is taken up in Chapter 5.

Section 4.6 extends stochastic calculus to multiple processes. Section $4.7$ discusses the Brownian bridge, which plays a useful role in Monte Carlo methods for pricing. We do not treat Monte Carlo methods in this text; we include the Brownian bridge only because it is a natural application of the stochastic calculus developed in the earlier sections.

\subsection{Itô's Integral for Simple Integrands}

We fix a positive number $T$ and seek to make sense of

$$
\int_{0}^{T} \Delta(t) d W(t) .
$$

The basic ingredients here are a Brownian motion $W(t), t \geq 0$, together with a filtration $\mathcal{F}(t), t \geq 0$, for this Brownian motion. We will let the integrand $\Delta(t)$ be an adapted stochastic process. Our reason for doing this is that $\Delta(t)$ will eventually be the position we take in an asset at time $t$, and this typically depends on the price path of the asset up to time $t$. Anything that depends on the path of a random process is itself random. Requiring $\Delta(t)$ to be adapted means that we require $\Delta(t)$ to be $\mathcal{F}(t)$-measurable for each $t \geq 0$. In other words, the information available at time $t$ is sufficient to evaluate $\Delta(t)$ at that time. When we are standing at time 0 and $t$ is strictly positive, $\Delta(t)$ is unknown to us. It is a random variable. When we get to time $t$, we have sufficient information to evaluate $\Delta(t)$; its randomness has been resolved.

Recall that increments of the Brownian motion after time $t$ are independent of $\mathcal{F}(t)$, and since $\Delta(t)$ is $\mathcal{F}(t)$-measurable, it must also be independent of these future Brownian increments. Positions we take in assets may depend on the price history of those assets, but they must be independent of the future increments of the Brownian motion that drives those prices.

The problem we face when trying to assign meaning to the Itô integral (4.2.1) is that Brownian motion paths cannot be differentiated with respect to time. If $g(t)$ is a differentiable function, then we can define

$$
\int_{0}^{T} \Delta(t) d g(t)=\int_{0}^{T} \Delta(t) g^{\prime}(t) d t,
$$

where the right-hand side is an ordinary (Lebesgue) integral with respect to time. This will not work for Brownian motion.

\subsubsection{Construction of the Integral}

To define the integral (4.2.1), Itô devised the following way around the nondifferentiability of the Brownian paths. We first define the Itô integral for simple integrands $\Delta(t)$ and then extend it to nonsimple integrands as a limit of the integral of simple integrands. We describe this procedure.

Let $\Pi=\left\{t_{0}, t_{1}, \ldots, t_{n}\right\}$ be a partition of $[0, T] ;$ i.e.,

$$
0=t_{0} \leq t_{1} \leq \cdots \leq t_{n}=T .
$$

Assume that $\Delta(t)$ is constant in $t$ on each subinterval $\left[t_{j}, t_{j+1}\right)$. Such a process $\Delta(t)$ is a simple process.

Figure 4.2.1 shows a single path of a simple process $\Delta(t)$. We shall always choose these simple processes, as shown in this figure, to take a value at a partition time $t_{j}$ and then hold it up to but not including the next partition time $t_{j+1}$. Although it is not apparent from Figure 4.2.1, the path shown depends on the same $\omega$ on which the path of the Brownian motion $W(t)$ (not shown) depends. If one were to choose a different $\omega$, there would be a different path of the Brownian motion and possibly a different path of $\Delta(t)$. However, the value of $\Delta(t)$ can depend only on the information available at time $t$. Since there is no information at time 0 , the value of $\Delta(0)$ must be the same for all paths, and hence the first piece of $\Delta(t)$, for $0 \leq t<t_{1}$, does not really depend on $\omega$. The value of $\Delta(t)$ on the second interval, $\left[t_{1}, t_{2}\right)$, can depend on observations made during the first time interval $\left[0, t_{1}\right)$. 

![](https://cdn.mathpix.com/cropped/2023_02_09_f310c50a6371948e10fbg-147.jpg?height=762&width=1098&top_left_y=163&top_left_x=96)

Fig. 4.2.1. A path of a simple process.

We shall think of the interplay between the simple process $\Delta(t)$ and the Brownian motion $W(t)$ in (4.2.1) in the following way. Regard $W(t)$ as the price per share of an asset at time $t$. (Since Brownian motion can take negative as well as positive values, it is not a good model of the price of a limitedliability asset such as a stock. For the sake of this illustration, we ignore that issue.) Think of $t_{0}, t_{1}, \ldots, t_{n-1}$ as the trading dates in the asset, and think of $\Delta\left(t_{0}\right), \Delta\left(t_{1}\right), \ldots, \Delta\left(t_{n-1}\right)$ as the position (number of shares) taken in the asset at each trading date and held to the next trading date. The gain from trading at each time $t$ is given by

$$
\begin{aligned}
& I(t)=\Delta\left(t_{0}\right)\left[W(t)-W\left(t_{0}\right)\right]=\Delta(0) W(t), \quad 0 \leq t \leq t_{1} \\
& I(t)=\Delta(0) W\left(t_{1}\right)+\Delta\left(t_{1}\right)\left[W(t)-W\left(t_{1}\right)\right], \quad t_{1} \leq t \leq t_{2} \\
& I(t)=\Delta(0) W\left(t_{1}\right)+\Delta\left(t_{1}\right)\left[W\left(t_{2}\right)-W\left(t_{1}\right)\right]+\Delta\left(t_{2}\right)\left[W(t)-W\left(t_{2}\right)\right] \\
& \quad t_{2} \leq t \leq t_{3}
\end{aligned}
$$

and so on. In general, if $t_{k} \leq t \leq t_{k+1}$, then

$$
I(t)=\sum_{j=0}^{k-1} \Delta\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]+\Delta\left(t_{k}\right)\left[W(t)-W\left(t_{k}\right)\right]
$$

The process $I(t)$ in (4.2.2) is the Itô integral of the simple process $\Delta(t)$, a fact that we write as 

$$
I(t)=\int_{0}^{t} \Delta(u) d W(u) .
$$

In particular, we can take $t=t_{n}=T$, and (4.2.2) provides a definition for the Itô integral (4.2.1). We have managed to define this integral not only for the upper limit of integration $T$ but also for every upper limit of integration $t$ between 0 and $T$.

\subsubsection{Properties of the Integral}

The Itô integral (4.2.2) is defined as the gain from trading in the martingale $W(t)$. A martingale has no tendency to rise or fall, and hence it is to be expected that $I(t)$, thought of as a process in its upper limit of integration $t$, also has no tendency to rise or fall. We formalize this observation by the next theorem and proof.

Theorem 4.2.1. The Itô integral defined by (4.2.2) is a martingale.

Proof: Let $0 \leq s \leq t \leq T$ be given. We shall assume that $s$ and $t$ are in different subintervals of the partition $\Pi$ (i.e., there are partition points $t_{\ell}$ and $t_{k}$ such that $t_{\ell}<t_{k}, s \in\left[t_{\ell}, t_{\ell+1}\right)$, and $\left.t \in\left[t_{k}, t_{k+1}\right)\right)$. If $s$ and $t$ are in the same subinterval, the following proof simplifies. Equation (4.2.2) may be rewritten as

$$
\begin{aligned}
I(t)= & \sum_{j=0}^{\ell-1} \Delta\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]+\Delta\left(t_{\ell}\right)\left[W\left(t_{\ell+1}\right)-W\left(t_{\ell}\right)\right] \\
& +\sum_{j=\ell+1}^{k-1} \Delta\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]+\Delta\left(t_{k}\right)\left[W(t)-W\left(t_{k}\right)\right] .
\end{aligned}
$$

We must show that $\mathbb{E}[I(t) \mid \mathcal{F}(s)]=I(s)$. We take the conditional expectation of each of the four terms on the right-hand side of (4.2.3). Every random variable in the first sum $\sum_{j=0}^{\ell-1} \Delta\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]$ is $\mathcal{F}(s)$-measurable because the latest time appearing in this sum is $t_{\ell}$ and $t_{\ell} \leq s$. Therefore,

$$
\mathbb{E}\left[\sum_{j=0}^{\ell-1} \Delta\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right] \mid \mathcal{F}(s)\right]=\sum_{j=0}^{\ell-1} \Delta\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]
$$

For the second term on the right-hand side of (4.2.3), we "take out what is known" (Theorem 2.3.2(ii)) and use the martingale property of $W$ to write

$$
\begin{aligned}
\mathbb{E}\left[\Delta\left(t_{\ell}\right)\left(W\left(t_{\ell+1}\right)-W\left(t_{\ell}\right)\right) \mid \mathcal{F}(s)\right] & =\Delta\left(t_{\ell}\right)\left(\mathbb{E}\left[W\left(t_{\ell+1}\right) \mid \mathcal{F}(s)\right]-W\left(t_{\ell}\right)\right) \\
& =\Delta\left(t_{\ell}\right)\left(W(s)-W\left(t_{\ell}\right)\right) .
\end{aligned}
$$

Adding (4.2.4) and (4.2.5), we obtain $I(s)$. It remains to show that the conditional expectations of the third and fourth terms on the right-hand side of (4.2.3) are zero. We will then have $\mathbb{E}[I(t) \mid \mathcal{F}(s)]=I(s)$

The summands in the third term are of the form $\Delta\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]$, where $t_{j} \geq t_{\ell+1}>s$. This permits us to use the following iterated conditioning trick, which is based on properties (iii) (iterated conditioning) and (ii) (taking out what is known) of Theorem 2.3.2:

$$
\begin{aligned}
\mathbb{E}\left\{\Delta\left(t_{j}\right)\right. & \left.\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right) \mid \mathcal{F}(s)\right\} \\
& =\mathbb{E}\left\{\mathbb{E}\left[\Delta\left(t_{j}\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right) \mid \mathcal{F}\left(t_{j}\right)\right] \mid \mathcal{F}(s)\right\} \\
& =\mathbb{E}\left\{\Delta\left(t_{j}\right)\left(\mathbb{E}\left[W\left(t_{j+1}\right) \mid \mathcal{F}\left(t_{j}\right)\right]-W\left(t_{j}\right)\right) \mid \mathcal{F}(s)\right\} \\
& =\mathbb{E}\left\{\Delta\left(t_{j}\right)\left(W\left(t_{j}\right)-W\left(t_{j}\right)\right) \mid \mathcal{F}(s)\right\}=0
\end{aligned}
$$

At the end, we have used the fact that $W$ is a martingale. Because the conditional expectation of each of the summands in the third term on the right-hand side of (4.2.3) is zero, the conditional expectation of the whole term is zero:

$$
\mathbb{E}\left\{\sum_{j=\ell+1}^{k-1} \Delta\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right] \mid \mathcal{F}(s)\right\}=0
$$

The fourth term on the right-hand side of (4.2.3) is treated like the summands in the third term, with the result that

$$
\begin{aligned}
\mathbb{E}\left\{\Delta\left(t_{k}\right)\right. & \left.\left(W(t)-W\left(t_{k}\right)\right) \mid \mathcal{F}(s)\right\} \\
& =\mathbb{E}\left\{\mathbb{E}\left[\Delta\left(t_{k}\right)\left(W(t)-W\left(t_{k}\right)\right) \mid \mathcal{F}\left(t_{k}\right)\right] \mid \mathcal{F}(s)\right\} \\
& =\mathbb{E}\left\{\Delta\left(t_{k}\right)\left(\mathbb{E}\left[W(t) \mid \mathcal{F}\left(t_{k}\right)\right]-W\left(t_{k}\right)\right) \mid \mathcal{F}(s)\right\} \\
& =\mathbb{E}\left\{\Delta\left(t_{k}\right)\left(W\left(t_{k}\right)-W\left(t_{k}\right)\right) \mid \mathcal{F}(s)\right\}=0 .
\end{aligned}
$$

This concludes the proof.

Because $I(t)$ is a martingale and $I(0)=0$, we have $\mathbb{E} I(t)=0$ for all $t \geq 0$. It follows that $\operatorname{Var} I(t)=\mathbb{E} I^{2}(t)$, a quantity that can be evaluated by the formula in the next theorem.

Theorem 4.2.2 (Itô isometry). The Itô integral defined by (4.2.2) satisfies

$$
\mathbb{E} I^{2}(t)=\mathbb{E} \int_{0}^{t} \Delta^{2}(u) d u
$$

Proof: To simplify the notation, we set $D_{j}=W\left(t_{j+1}\right)-W\left(t_{j}\right)$ for $j=$ $0, \ldots, k-1$ and $D_{k}=W(t)-W\left(t_{k}\right)$ so that (4.2.2) may be written as $I(t)=$ $\sum_{j=0}^{k} \Delta\left(t_{j}\right) D_{j}$ and

$$
I^{2}(t)=\sum_{j=0}^{k} \Delta^{2}\left(t_{j}\right) D_{j}^{2}+2 \sum_{0 \leq i<j \leq k} \Delta\left(t_{i}\right) \Delta\left(t_{j}\right) D_{i} D_{j}
$$

We first show that the expected value of each of the cross terms is zero. For $i<$ $j$, the random variable $\Delta\left(t_{i}\right) \Delta\left(t_{j}\right) D_{i}$ is $\mathcal{F}\left(t_{j}\right)$-measurable, while the Brownian increment $D_{j}$ is independent of $\mathcal{F}\left(t_{j}\right)$. Furthermore, $\mathbb{E} D_{j}=0$. Therefore,

$$
\mathbb{E}\left[\Delta\left(t_{i}\right) \Delta\left(t_{j}\right) D_{i} D_{j}\right]=\mathbb{E}\left[\Delta\left(t_{i}\right) \Delta\left(t_{j}\right) D_{i}\right] \cdot \mathbb{E} D_{j}=\mathbb{E}\left[\Delta\left(t_{i}\right) \Delta\left(t_{j}\right) D_{i}\right] \cdot 0=0
$$

We next consider the square terms $\Delta^{2}\left(t_{j}\right) D_{j}^{2}$. The random variable $\Delta^{2}\left(t_{j}\right)$ is $\mathcal{F}\left(t_{j}\right)$-measurable, and the squared Brownian increment $D_{j}^{2}$ is independent of $\mathcal{F}\left(t_{j}\right)$. Furthermore, $\mathbb{E} D_{j}^{2}=t_{j+1}-t_{j}$ for $j=0, \ldots, k-1$ and $\mathbb{E} D_{k}^{2}=t-t_{k}$. Therefore,

$$
\begin{aligned}
\mathbb{E} I^{2}(t) & =\sum_{j=0}^{k} \mathbb{E}\left[\Delta^{2}\left(t_{j}\right) D_{j}^{2}\right]=\sum_{j=1}^{k} \mathbb{E} \Delta^{2}\left(t_{j}\right) \cdot \mathbb{E} D_{j}^{2} \\
& =\sum_{j=0}^{k-1} \mathbb{E} \Delta^{2}\left(t_{j}\right)\left(t_{j+1}-t_{j}\right)+\mathbb{E} \Delta^{2}\left(t_{k}\right)\left(t-t_{k}\right) .
\end{aligned}
$$

But $\Delta\left(t_{j}\right)$ is constant on the interval $\left[t_{j}, t_{j+1}\right)$, and hence $\Delta^{2}\left(t_{j}\right)\left(t_{j+1}-t_{j}\right)=$ $\int_{t_{j}}^{t_{j+1}} \Delta^{2}(u) d u$. Similarly, $\Delta^{2}\left(t_{k}\right)\left(t-t_{k}\right)=\int_{t_{k}}^{t} \Delta^{2}(u) d u$. We may thus continue $(4.2 .7)$ to obtain

$$
\begin{aligned}
\mathbb{E} I^{2}(t) & =\sum_{j=0}^{k-1} \mathbb{E} \int_{t_{j}}^{t_{j+1}} \Delta^{2}(u) d u+\mathbb{E} \int_{t_{k}}^{t} \Delta^{2}(u) d u \\
& =\mathbb{E}\left[\sum_{j=0}^{k-1} \int_{t_{j}}^{t_{j+1}} \Delta^{2}(u) d u+\int_{t_{k}}^{t} \Delta^{2}(u) d u\right]=\mathbb{E} \int_{0}^{t} \Delta^{2}(u) d u .
\end{aligned}
$$

Finally, we turn to the quadratic variation of the Itô integral $I(t)$ thought of as a process in its upper limit of integration $t$. Brownian motion accumulates quadratic variation at rate one per unit time. However, Brownian motion is scaled in a time- and path-dependent way by the integrand $\Delta(u)$ as it enters the Itô integral $I(t)=\int_{0}^{t} \Delta(u) d B(u)$. Because increments are squared in the computation of quadratic variation, the quadratic variation of Brownian motion will be scaled by $\Delta^{2}(u)$ as it enters the Itô integral. The following theorem gives the precise statement. Theorem 4.2.3. The quadratic variation accumulated up to time $t$ by the Itô integral (4.2.2) is

$$
[I, I](t)=\int_{0}^{t} \Delta^{2}(u) d u
$$

Proof: We first compute the quadratic variation accumulated by the Itô integral on one of the subintervals $\left[t_{j}, t_{j+1}\right]$ on which $\Delta(u)$ is constant. For this, we choose partition points

$$
t_{j}=s_{0}<s_{1}<\cdots<s_{m}=t_{j+1}
$$

and consider

$$
\begin{aligned}
\sum_{i=0}^{m-1}\left[I\left(s_{i+1}\right)-I\left(s_{i}\right)\right]^{2} & =\sum_{i=0}^{m-1}\left[\Delta\left(t_{j}\right)\left(W\left(s_{i+1}\right)-W\left(s_{i}\right)\right)\right]^{2} \\
& =\Delta^{2}\left(t_{j}\right) \sum_{i=0}^{m-1}\left(W\left(s_{i+1}\right)-W\left(s_{i}\right)\right)^{2}
\end{aligned}
$$

As $m \rightarrow \infty$ and the step size $\max _{i=0, \ldots, m-1}\left(s_{i+1}-s_{i}\right)$ approaches zero, the term $\sum_{i=0}^{m-1}\left(W\left(s_{i+1}\right)-W\left(s_{i}\right)\right)^{2}$ converges to the quadratic variation accumulated by Brownian motion between times $t_{j}$ and $t_{j+1}$, which is $t_{j+1}-t_{j}$. Therefore, the limit of $(4.2 .9)$, which is the quadratic variation accumulated by the Itô integral between times $t_{j}$ and $t_{j+1}$, is

$$
\Delta^{2}\left(t_{j}\right)\left(t_{j+1}-t_{j}\right)=\int_{t_{j}}^{t_{j+1}} \Delta^{2}(u) d u,
$$

where again we have used the fact that $\Delta(u)$ is constant for $t_{j} \leq u<t_{j+1}$. Analogously, the quadratic variation accumulated by the Itô integral between times $t_{k}$ and $t$ is $\int_{t_{k}}^{t} \Delta^{2}(u) d u$. Adding up all these pieces, we obtain (4.2.8).

In Theorems $4.2 .2$ and $4.2 .3$, we finally see how the quadratic variation and the variance of a process can differ. The quadratic variation is computed path-by-path, and the result can depend on the path. If along one path of the Brownian motion we choose large positions $\Delta(u)$, the Itô integral will have a large quadratic variation. Along a different path, we could choose small positions $\Delta(u)$ and the Itô integral would have a small quadratic variation. The quadratic variation can be regarded as a measure of risk, and it depends on the size of the positions we take. The variance of $I(t)$ is an average over all possible paths of the quadratic variation. Because it is the expectation of something, it cannot be random. As an average over all possible paths, realized and unrealized, it is a more theoretical concept than quadratic variation. We emphasize here that what we are calling variance is not the empirical variance. Empirical (or sample) variance is computed from a realized path and is an estimator of the theoretical variance we are discussing. The empirical variance is sometimes carelessly called variance, which creates the possibility of confusion.

Finally, we recall the equation (3.4.10), $d W(t) d W(t)=d t$, of Remark 3.4.4. We interpret this equation as the statement that Brownian motion accumulates quadratic variation at rate one per unit time. It is another way of writing $[W, W](t)=t, t \geq 0$. The Itô integral formula $I(t)=\int_{0}^{t} \Delta(u) d W(u)$ can be written in differential form as $d I(t)=\Delta(t) d W(t)$, and we can then use (3.4.10) to square $d I(t)$ :

$$
d I(t) d I(t)=\Delta^{2}(t) d W(t) d W(t)=\Delta^{2}(t) d t .
$$

This equation says that the Itô integral $I(t)$ accumulates quadratic variation at rate $\Delta^{2}(t)$ per unit time. The rate of accumulation is typically both timeand path-dependent. Equation (4.2.10) is another way of reporting the result of Theorem 4.2.3.

Remark 4.2.4 (on notation). The notations

$$
I(t)=\int_{0}^{t} \Delta(u) d W(u)
$$

and

$$
d I(t)=\Delta(t) d W(t)
$$

mean almost the same thing, although the second is probably more intuitive. Equation (4.2.11) has the precise meaning given by (4.2.2). Equation (4.2.12) has the imprecise meaning that when we move forward a little bit in time from time $t$, the change in the Itô integral $I$ is $\Delta(t)$ times the change in the Brownian motion $W$. It also has a precise meaning, which one obtains by integrating both sides, remembering to put in a constant of integration $I(0)$ :

$$
I(t)=I(0)+\int_{0}^{t} \Delta(u) d W(u) .
$$

We say that (4.2.12) is the differential form of (4.2.13) and that (4.2.13) is the integral form of (4.2.12). These two equations mean exactly the same thing.

The only difference between (4.2.11) and (4.2.13), and hence the only difference between (4.2.11) and (4.2.12), is that (4.2.11) specifies the initial condition $I(0)=0$, whereas (4.2.12) and (4.2.13) permit $I(0)$ to be any arbitrary constant.

\subsection{Itô's Integral for General Integrands}

In this section, we define the Itô integral $\int_{0}^{T} \Delta(t) d W(t)$ for integrands $\Delta(t)$ that are allowed to vary continuously with time and also to jump. In particular, we no longer assume that $\Delta(t)$ is a simple process as shown in Figure 4.2.1. We do assume that $\Delta(t), t \geq 0$, is adapted to the filtration $\mathcal{F}(t), t \geq 0$. We also assume the square-integrability condition

$$
\mathbb{E} \int_{0}^{T} \Delta^{2}(t) d t<\infty .
$$

In order to define $\int_{0}^{T} \Delta(t) d W(t)$, we approximate $\Delta(t)$ by simple processes. Figure 4.3.1 suggests how this can be done. In that figure, the continuously varying $\Delta(t)$ is shown as a solid line and the approximating simple integrand is dashed. Notice that $\Delta(t)$ is allowed to jump. The approximating simple integrand is constructed by choosing a partition $0=t_{0}<t_{1}<t_{2}<t_{3}<t_{4}$, setting the approximating simple process equal to $\Delta\left(t_{j}\right)$ at each $t_{j}$, and then holding the simple process constant over the subinterval $\left[t_{j}, t_{j+1}\right)$. As the maximal step size of the partition approaches zero, the approximating integrand will become a better and better approximation of the continuously varying one.

![](https://cdn.mathpix.com/cropped/2023_02_09_f310c50a6371948e10fbg-153.jpg?height=772&width=1102&top_left_y=875&top_left_x=90)

Fig. 4.3.1. Approximating a continuously varying integrand.

In general, then, it is possible to choose a sequence $\Delta_{n}(t)$ of simple processes such that as $n \rightarrow \infty$ these processes converge to the continuously varying $\Delta(t)$. By "converge," we mean that

$$
\lim _{n \rightarrow \infty} \mathbb{E} \int_{0}^{T}\left|\Delta_{n}(t)-\Delta(t)\right|^{2} d t=0
$$

For each $\Delta_{n}(t)$, the Itô integral $\int_{0}^{t} \Delta_{n}(u) d W(u)$ has already been defined for $0 \leq t \leq T$. We define the Itô integral for the continuously varying integrand $\Delta(t)$ by the formula 1

$$
\int_{0}^{t} \Delta(u) d W(u)=\lim _{n \rightarrow \infty} \int_{0}^{t} \Delta_{n}(u) d W(u), \quad 0 \leq t \leq T .
$$

This integral inherits the properties of Itô integrals of simple processes. We summarize these in the next theorem.

Theorem 4.3.1. Let $T$ be a positive constant and let $\Delta(t), 0 \leq t \leq T$, be an adapted stochastic process that satisfies (4.3.1). Then $I(t)=\int_{0}^{t} \Delta(u) d W(u)$ defined by (4.3.3) has the following properties.

(i) (Continuity) As a function of the upper limit of integration the paths of $I(t)$ are continuous.

(ii) (Adaptivity) For each $t, I(t)$ is $\mathcal{F}(t)$-measurable.

(üi) (Linearity) If $I(t)=\int_{0}^{t} \Delta(u) d W(u)$ and $J(t)=\int_{0}^{t} \Gamma(u) d W(u)$, then $I(t) \pm J(t)=\int_{0}^{t}(\Delta(u) \pm \Gamma(u)) d W(u)$; furthermore, for every constant $c$, $c I(t)=\int_{0}^{t} c \Delta(u) d W(u)$.

(iv) (Martingale) $I(t)$ is a martingale.

(v) (Itô isometry) $\mathbb{E} I^{2}(t)=\mathbb{E} \int_{0}^{t} \Delta^{2}(u) d u$.

(vi) (Quadratic variation) $[I, I](t)=\int_{0}^{t} \Delta^{2}(u) d u$.

Example 4.3.2. We compute $\int_{0}^{T} W(t) d W(t)$. To do that, we choose a large integer $n$ and approximate the integrand $\Delta(t)=W(t)$ by the simple process

$$
\Delta_{n}(t)= \begin{cases}W(0)=0 & \text { if } 0 \leq t<\frac{T}{n}, \\ W\left(\frac{T}{n}\right) & \text { if } \frac{T}{n} \leq t<\frac{2 T}{n}, \\ \vdots \\ W\left(\frac{(n-1) T}{n}\right) & \text { if } \frac{(n-1) T}{n} \leq t<T,\end{cases}
$$

as shown in Figure 4.3.2. Then $\lim _{n \rightarrow \infty} \mathbb{E} \int_{0}^{T}\left|\Delta_{n}(t)-W(t)\right|^{2} d t=0$. By definition,

$$
\begin{aligned}
& \int_{0}^{T} W(t) d W(t)=\lim _{n \rightarrow \infty} \int_{0}^{T} \Delta_{n}(t) d W(t) \\
& \quad=\lim _{n \rightarrow \infty} \sum_{j=0}^{n-1} W\left(\frac{j T}{n}\right)\left[W\left(\frac{(j+1) T}{n}\right)-W\left(\frac{j T}{n}\right)\right] .
\end{aligned}
$$

${ }^{1}$ For each $t$, the limit in (4.3.3) exists because $I_{n}(t)=\int_{0}^{t} \Delta_{n}(u) d W(u)$ is a Cauchy sequence in $L_{2}(\Omega, \mathcal{F}, \mathbb{P})$. This is because of Itô's isometry (Theorem 4.2.2), which yields $\mathbb{E}\left(I_{n}(t)-I_{m}(t)\right)^{2}=\mathbb{E} \int_{0}^{t}\left|\Delta_{n}(u)-\Delta_{m}(u)\right|^{2} d u$. As a consequence of (4.3.2), the right-hand side has limit zero as $n$ and $m$ approach infinity. 

![](https://cdn.mathpix.com/cropped/2023_02_09_f310c50a6371948e10fbg-155.jpg?height=985&width=1066&top_left_y=100&top_left_x=160)

Fig. 4.3.2. Simple process approximating Brownian motion.

To simplify notation, we denote $W_{j}=W\left(\frac{j T}{n}\right)$. As a precursor to evaluating the limit in (4.3.4), we work out equation (4.3.5) below. The second equality in (4.3.5) is obtained by making the change of index $k=j+1$ in the first sum. The third equality uses the fact that $W_{0}=W(0)=0$. We have

$$
\begin{aligned}
\frac{1}{2} \sum_{j=0}^{n-1}\left(W_{j+1}-W_{j}\right)^{2} & =\frac{1}{2} \sum_{j=0}^{n-1} W_{j+1}^{2}-\sum_{j=0}^{n-1} W_{j} W_{j+1}+\frac{1}{2} \sum_{j=0}^{n-1} W_{j}^{2} \\
& =\frac{1}{2} \sum_{k=1}^{n} W_{k}^{2}-\sum_{j=0}^{n-1} W_{j} W_{j+1}+\frac{1}{2} \sum_{j=0}^{n-1} W_{j}^{2} \\
& =\frac{1}{2} W_{n}^{2}+\frac{1}{2} \sum_{k=0}^{n-1} W_{k}^{2}-\sum_{j=0}^{n-1} W_{j} W_{j+1}+\frac{1}{2} \sum_{j=0}^{n-1} W_{j}^{2} \\
& =\frac{1}{2} W_{n}^{2}+\sum_{j=0}^{n-1} W_{j}^{2}-\sum_{j=0}^{n-1} W_{j} W_{j+1} \\
& =\frac{1}{2} W_{n}^{2}+\sum_{j=0}^{n-1} W_{j}\left(W_{j}-W_{j+1}\right)
\end{aligned}
$$

From (4.3.5), we conclude that

$$
\sum_{j=0}^{n-1} W_{j}\left(W_{j+1}-W_{j}\right)=\frac{1}{2} W_{n}^{2}-\frac{1}{2} \sum_{j=0}^{n-1}\left(W_{j+1}-W_{j}\right)^{2} .
$$

In the original notation, this is

$$
\begin{aligned}
\sum_{j=0}^{n-1} W & \left(\frac{j T}{n}\right)\left[W\left(\frac{(j+1) T}{n}\right)-W\left(\frac{j T}{n}\right)\right] \\
& =\frac{1}{2} W^{2}(T)-\frac{1}{2} \sum_{j=0}^{n-1}\left[W\left(\frac{(j+1) T}{n}\right)-W\left(\frac{j T}{n}\right)\right]^{2}
\end{aligned}
$$

Letting $n \rightarrow \infty$ in (4.3.4) and using this equation, we get

$$
\int_{0}^{T} W(t) d W(t)=\frac{1}{2} W^{2}(T)-\frac{1}{2}[W, W](T)=\frac{1}{2} W^{2}(T)-\frac{1}{2} T .
$$

We contrast (4.3.6) with ordinary calculus. If $g$ is a differentiable function with $g(0)=0$, then

$$
\int_{0}^{T} g(t) d g(t)=\int_{0}^{T} g(t) g^{\prime}(t) d t=\left.\frac{1}{2} g^{2}(t)\right|_{0} ^{T}=\frac{1}{2} g^{2}(T)
$$

The extra term $-\frac{1}{2} T$ in (4.3.6) comes from the nonzero quadratic variation of Brownian motion and the way we constructed the Itô integral, always evaluating the integrand at the left-hand endpoint of the subinterval (see the right-hand side of (4.3.4)). If we were instead to evaluate at the midpoint, replacing the right-hand side of $(4.3 .4)$ by

$$
\lim _{n \rightarrow \infty} \sum_{j=0}^{n-1} W\left(\frac{\left(j+\frac{1}{2}\right) T}{n}\right)\left[W\left(\frac{(j+1) T}{n}\right)-W\left(\frac{j T}{n}\right)\right]
$$

then we would not have gotten this term (see Exercise 4.4). The integral obtained by making this replacement is called the Stratonovich integral, and the ordinary rules of calculus apply to it. However, it is inappropriate for finance. In finance, the integrand represents a position in an asset and the integrator represents the price of that asset. We cannot decide at 1:00 p.m. which position we took at 9:00 a.m. We must decide the position at the beginning of each time interval, and the Itô integral is the limit of the gain achieved by that kind of trading as the time between trades approaches zero.

For functions $g(t)$ that have a derivative, integrals such as $\int_{0}^{t} g(t) d g(t)$ are not sensitive to this distinction (i.e., the Itô integral and Stratonovich integral approximations have the same limit, which is $\frac{1}{2} g^{2}(T)$ ). For functions that have a nonzero quadratic variation, integrals are sensitive to where in the subintervals the approximating integrands are evaluated. The upper limit of integration $T$ in (4.3.6) is arbitrary and can be replaced by any $t \geq 0$. In other words,

$$
\int_{0}^{t} W(u) d W(u)=\frac{1}{2} W^{2}(t)-\frac{1}{2} t, \quad t \geq 0 .
$$

Theorem 4.3.1(iv) guarantees that $\int_{0}^{t} W(u) d W(u)$ is a martingale and hence has constant expectation. At $t=0$, this martingale is 0 , and hence its expectation must always be zero. This is indeed the case because $\mathbb{E} W^{2}(t)=t$. If the term $-\frac{1}{2} t$ were not present, we would not have a martingale.

\subsection{Itô-Doeblin Formula}

The addition of Doeblin's name to what has traditionally been called the Itô formula is explained in the Notes, Section 4.9.

\subsubsection{Formula for Brownian Motion}

We want a rule to "differentiate" expressions of the form $f(W(t))$, where $f(x)$ is a differentiable function and $W(t)$ is a Brownian motion. If $W(t)$ were also differentiable, then the chain rule from ordinary calculus would give

$$
\frac{d}{d t} f(W(t))=f^{\prime}(W(t)) W^{\prime}(t),
$$

which could be written in differential notation as

$$
d f(W(t))=f^{\prime}(W(t)) W^{\prime}(t) d t=f^{\prime}(W(t)) d W(t) .
$$

Because $W$ has nonzero quadratic variation, the correct formula has an extra term, namely,

$$
d f(W(t))=f^{\prime}(W(t)) d W(t)+\frac{1}{2} f^{\prime \prime}(W(t)) d t .
$$

This is the Itô-Doeblin formula in differential form. Integrating this, we obtain the Itô-Doeblin formula in integral form:

$$
f(W(t))-f(W(0))=\int_{0}^{t} f^{\prime}(W(u)) d W(u)+\frac{1}{2} \int_{0}^{t} f^{\prime \prime}(W(u)) d u .
$$

The mathematically meaningful form of the Itô-Doeblin formula is the integral form (4.4.2). This is because we have precise definitions for both terms appearing on the right-hand side. The first, $\int_{0}^{t} f^{\prime}(W(u)) d W(u)$, is an Itô integral, defined in the previous section. The second, $\int_{0}^{t} f^{\prime \prime}(W(u)) d u$, is an ordinary (Lebesgue) integral with respect to the time variable. For pencil and paper computations, the more convenient form of the ItôDoeblin formula is the differential form (4.4.1). There is an intuitive meaning but no precise definition for the terms $d f(W(t)), d W(t)$, and $d t$ appearing in this formula. The intuitive meaning is that $d f(W(t))$ is the change in $f(W(t))$ when $t$ changes a "little bit" $d t, d W(t)$ is the change in the Brownian motion when $t$ changes a "little bit" $d t$, and the whole formula is exact only when the "little bit" is "infinitesimally small." Because there is no precise definition for "little bit" and "infinitesimally small," we rely on (4.4.2) to give precise meaning to (4.4.1).

The relationship between (4.4.1) and (4.4.2) is similar to that developed in ordinary calculus to assist in changing variables in an integral. If asked to compute the indefinite integral $\int f(u) f^{\prime}(u) d u$, we might make the change of variable $v=f(u)$ and write $d v=f^{\prime}(u) d u$, so that the indefinite integral becomes $\int v d v$, which is $\frac{1}{2} v^{2}+C=\frac{1}{2} f^{2}(u)+C$, where $C$ is a constant of integration. The final formula

$$
\int f(u) f^{\prime}(u) d u=\frac{1}{2} f^{2}(u)+C
$$

is correct, as can be verified by differentiating $\frac{1}{2} f^{2}(u)+C$ to get $f(u) f^{\prime}(u)$. We do not attempt to give precise definitions to the terms $d v$ and $d u$ appearing in the equation $d v=f^{\prime}(u) d u$ used in deriving it.

We formalize the preceding discussion with a theorem that provides a formula slightly more general than (4.4.2) in that it allows $f$ to be a function of both $t$ and $x$.

Theorem 4.4.1 (Itô-Doeblin formula for Brownian motion). Let $f(t, x)$ be a function for which the partial derivatives $f_{t}(t, x), f_{x}(t, x)$, and $f_{x x}(t, x)$ are defined and continuous, and let $W(t)$ be a Brownian motion. Then, for every $T \geq 0$,

$$
\begin{aligned}
f(T, W(T))=f & (0, W(0))+\int_{0}^{T} f_{t}(t, W(t)) d t \\
& +\int_{0}^{T} f_{x}(t, W(t)) d W(t)+\frac{1}{2} \int_{0}^{T} f_{x x}(t, W(t)) d t
\end{aligned}
$$

SKETCH OF PROOF: We first show why (4.4.3) holds when $f(x)=\frac{1}{2} x^{2}$. In this case, $f^{\prime}(x)=x$ and $f^{\prime \prime}(x)=1$. Let $x_{j+1}$ and $x_{j}$ be numbers. Taylor's formula implies

$$
f\left(x_{j+1}\right)-f\left(x_{j}\right)=f^{\prime}\left(x_{j}\right)\left(x_{j+1}-x_{j}\right)+\frac{1}{2} f^{\prime \prime}\left(x_{j}\right)\left(x_{j+1}-x_{j}\right)^{2} .
$$

In this case, Taylor's formula to second order is exact (there is no remainder term) because $f^{\prime \prime \prime}$ and all higher derivatives of $f$ are zero. We return to this matter later. Fix $T>0$, and let $\Pi=\left\{t_{0}, t_{1}, \ldots, t_{n}\right\}$ be a partition of $[0, T]$ (i.e., $0=$ $\left.t_{0}<t_{1}<\cdots<t_{n}=T\right)$. We are interested in the difference between $f(W(0))$ and $f(W(T))$. This change in $f(W(t))$ between times $t=0$ and $t=T$ can be written as the sum of the changes in $f(W(t))$ over each of the subintervals $\left[t_{j}, t_{j+1}\right]$. We do this and then use Taylor's formula (4.4.4) with $x_{j}=W\left(t_{j}\right)$ and $x_{j+1}=W\left(t_{j+1}\right)$ to obtain

$$
\begin{aligned}
f(W(T))-f(W(0))= & \sum_{j=0}^{n-1}\left[f\left(W\left(t_{j+1}\right)\right)-f\left(W\left(t_{j}\right)\right)\right] \\
= & \sum_{j=0}^{n-1} f^{\prime}\left(W\left(t_{j}\right)\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right] \\
& +\frac{1}{2} \sum_{j=0}^{n-1} f^{\prime \prime}\left(W\left(t_{j}\right)\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]^{2} .
\end{aligned}
$$

For the function $f(x)=\frac{1}{2} x^{2}$, the right-hand side of (4.4.5) is

$$
\sum_{j=0}^{n-1} W\left(t_{j}\right)\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]+\frac{1}{2} \sum_{j=0}^{n-1}\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]^{2} .
$$

If we let $\|\Pi\| \rightarrow 0$, the left-hand side of (4.4.5) is unaffected and the terms on the right-hand side converge to an Itô integral and one-half of the quadratic variation of Brownian motion, respectively:

$$
\begin{aligned}
& f(W(T))-f(W(0)) \\
& =\lim _{\|\Pi\| \rightarrow 0} \sum_{j=0}^{n-1} W\left(t_{j}\right)\left[W\left(v_{j+1}\right)-W\left(t_{j}\right)\right]+\lim _{\|\Pi\| \rightarrow 0} \frac{1}{2} \sum_{j=0}^{n-1}\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]^{2} \\
& =\int_{0}^{T} W(t) d W(t)+\frac{1}{2} T \\
& =\int_{0}^{T} f^{\prime}(W(t)) d W(t)+\frac{1}{2} \int_{0}^{T} f^{\prime \prime}(W(t)) d t .
\end{aligned}
$$

This is the Itô-Doeblin formula in integral form for the function $f(x)=\frac{1}{2} x^{2}$.

If instead of the quadratic function $f(x)=\frac{1}{2} x^{2}$ we had a general function $f(x)$, then in (4.4.5) we would have also gotten a sum of terms containing $\left[W\left(t_{j+1}\right)-W\left(t_{j}\right)\right]^{3}$. But according to Exercise $3.4$ of Chapter 3, $\sum_{j=0}^{n-1}\left|W\left(t_{j+1}\right)-W\left(t_{j}\right)\right|^{3}$ has limit zero as $\|\Pi\| \rightarrow 0$. Therefore, this term would make no contribution to the final answer.

If we take a function $f(t, x)$ of both the time variable $t$ and the variable $x$, then Taylor's Theorem says that 

$$
\begin{aligned}
& f\left(t_{j+1}, x_{j+1}\right)-f\left(t_{j}, x_{j}\right) \\
& =f_{t}\left(t_{j}, x_{j}\right)\left(t_{j+1}-t_{j}\right)+f_{x}\left(t_{j}, x_{j}\right)\left(x_{j+1}-x_{j}\right) \\
& \quad+\frac{1}{2} f_{x x}\left(t_{j}, x_{j}\right)\left(x_{j+1}-x_{j}\right)^{2}+f_{t x}\left(t_{j}, x_{j}\right)\left(t_{j+1}-t_{j}\right)\left(x_{j+1}-x_{j}\right) \\
& \quad+\frac{1}{2} f_{t t}\left(t_{j}, x_{j}\right)\left(t_{j+1}-t_{j}\right)^{2}+\text { higher-order terms. }
\end{aligned}
$$

We replace $x_{j}$ by $W\left(t_{j}\right)$, replace $x_{j+1}$ by $W\left(t_{j+1}\right)$, and sum:

$$
\begin{aligned}
& f(T, W(T))-f(0, W(0)) \\
& =\sum_{j=0}^{n-1}\left[f\left(t_{j+1}, W\left(t_{j+1}\right)\right)-f\left(t_{j}, W\left(t_{j}\right)\right)\right] \\
& =\sum_{j=0}^{n-1} f_{t}\left(t_{j}, W\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)+\sum_{j=0}^{n-1} f_{x}\left(t_{j}, W\left(t_{j}\right)\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right) \\
& \quad+\frac{1}{2} \sum_{j=0}^{n-1} f_{x x}\left(t_{j}, W\left(t_{j}\right)\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right)^{2} \\
& \quad+\sum_{j=0}^{n-1} f_{t x}\left(t_{j}, W\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right) \\
& \quad+\frac{1}{2} \sum_{j=0}^{n-1} f_{t t}\left(t_{j}, W\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)^{2}+\text { higher-order terms. }
\end{aligned}
$$

When we take the limit as $\|\Pi\| \rightarrow 0$, the left-hand side of (4.4.9) is unaffected. The first term on the right-hand side of (4.4.9) contributes the ordinary (Lebesgue) integral

$$
\lim _{\|\Pi\| \rightarrow 0} \sum_{j=0}^{n-1} f_{t}\left(t_{j}, W\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)=\int_{0}^{T} f_{t}(t, W(t)) d t
$$

to the final answer. As $\|\Pi\| \rightarrow 0$, the second term contributes the Itô integral $\int_{0}^{T} f_{x}(t, W(t)) d W(t)$. The third term contributes another ordinary (Lebesgue) integral, $\frac{1}{2} \int_{0}^{T} f_{x x}(t, W(t)) d t$, similar to the way we obtained this integral in (4.4.7). In other words, in the third term we can replace $\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right)^{2}$ by $t_{j+1}-t_{j}$. This is not an exact substitution, but when we sum the terms this substitution gives the correct limit as $\|\Pi\| \rightarrow 0$. See Remark 3.4.4 for more discussion of this point. With this substitution, the third term on the right-hand side of (4.4.9) contributes $\frac{1}{2} \int_{0}^{T} f_{x x}(t, W(t)) d t$. These limits of the first three terms appear on the right-hand side of (4.4.3). The fourth and fifth terms contribute zero. Indeed, for the fourth term, we observe that 

$$
\begin{aligned}
& \lim _{\|\Pi\| \rightarrow 0}\left|\sum_{j=0}^{n-1} f_{t x}\left(t_{j}, W\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right)\right| \\
& \leq \lim _{\|\Pi\| \rightarrow 0} \sum_{j=0}^{n-1}\left|f_{t x}\left(t_{j}, W\left(t_{j}\right)\right)\right| \cdot\left(t_{j+1}-t_{j}\right) \cdot\left|W\left(t_{j+1}\right)-W\left(t_{j}\right)\right| \\
& \leq \lim _{\|\Pi\| \rightarrow 0} \max _{0 \leq k \leq n-1}\left|W\left(t_{k+1}\right)-W\left(t_{k}\right)\right| \cdot \lim _{\|\Pi\| \rightarrow 0} \sum_{j=0}^{n-1}\left|f_{t x}\left(t_{j}, W\left(t_{j}\right)\right)\right|\left(t_{j+1}-t_{j}\right) \\
& =0 \cdot \int_{0}^{T} \mid f_{t x}(t, W(t)) d t=0 .
\end{aligned}
$$

The fifth term is treated similarly:

$$
\begin{aligned}
& \lim _{\|\Pi\| \rightarrow 0}\left|\frac{1}{2} \sum_{j=0}^{n-1} f_{t t}\left(t_{j}, W\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)^{2}\right| \\
& \leq \lim _{\|\Pi\| \rightarrow 0} \frac{1}{2} \sum_{j=0}^{n-1}\left|f_{t t}\left(t_{j}, W\left(t_{j}\right)\right)\right| \cdot\left(t_{j+1}-t_{j}\right)^{2} \\
& \leq \frac{1}{2} \lim _{\|\Pi\| \rightarrow 0} \max _{0 \leq n-1}\left(t_{k+1}-t_{k}\right) \cdot \lim _{\|\Pi\| \rightarrow 0} \sum_{j=0}^{n-1}\left|f_{t t}\left(t_{j}, W\left(t_{j}\right)\right)\right|\left(t_{j+1}-t_{j}\right) \\
&=\frac{1}{2} \cdot 0 \cdot \int_{0}^{T} f_{t t}(t, W(t)) d t=0 .
\end{aligned}
$$

The higher-order terms likewise contribute zero to the final answer.

Remark 4.4.2. The fact that the sum (4.4.10) of terms containing the product $\left(t_{j+1}-t_{j}\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right)$ has limit zero can be informally recorded by the formula $d t d W(t)=0$. Similarly, the sum (4.4.11) of terms containing $\left(t_{j+1}-t_{j}\right)^{2}$ also has limit zero, and this can be recorded by the formula $d t d t=0$. We can write these terms if we like in the Itô-Doeblin formula, so that in differential form it becomes

$$
\begin{aligned}
& d f(t, W(t)) \\
& =f_{t}(t, W(t)) d t+f_{x}(t, W(t)) d W(t)+\frac{1}{2} f_{x x}(t, W(t)) d W(t) d W(t) \\
& \quad+f_{t x}(t, W(t)) d t d W(t)+\frac{1}{2} f_{x x}(t, W(t)) d t d t
\end{aligned}
$$

but

$$
d W(t) d W(t)=d t, \quad d t d W(t)=d W(t) d t=0, \quad d t d t=0,
$$

and the Itô-Doeblin formula in differential form simplifies to

$$
d f(t, W(t))=f_{t}(t, W(t)) d t+f_{x}(t, W(t)) d W(t)+\frac{1}{2} f_{x x}(t, W(t)) d t .
$$

In Figure 4.4.1, we illustrate the Taylor series approximation of the difference $f\left(W\left(t_{j+1}\right)\right)-f\left(W\left(t_{j}\right)\right)$ for a function $f(x)$ that does not depend on $t$. The first-order approximation, which is $f^{\prime}\left(W\left(t_{j}\right)\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right)$, has an error due to the convexity of the function $f(x)$. Most of this error is removed by adding in the second-order term $\frac{1}{2} f^{\prime \prime}\left(W\left(t_{j}\right)\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right)^{2}$, which captures the curvature of the function $f(x)$ at $x=W\left(t_{j}\right)$.

![](https://cdn.mathpix.com/cropped/2023_02_09_f310c50a6371948e10fbg-162.jpg?height=617&width=1070&top_left_y=439&top_left_x=170)

Fig. 4.4.1. Taylor approximation to $f\left(W\left(t_{j+1}\right)\right)-f\left(W\left(t_{j}\right)\right)$.

In other words,

$f\left(W\left(t_{j+1}\right)\right)-f\left(W\left(t_{j}\right)\right)=f^{\prime}\left(W\left(t_{j}\right)\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right)+$ small error,

and

$$
\begin{aligned}
& f\left(W\left(t_{j+1}\right)\right)-f\left(W\left(t_{j}\right)\right)=f^{\prime}\left(W\left(t_{j}\right)\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right) \\
&+\frac{1}{2} f^{\prime \prime}\left(W\left(t_{j}\right)\right)\left(W\left(t_{j+1}\right)-W\left(t_{j}\right)\right)^{2} \\
&+\text { smaller error. }
\end{aligned}
$$

In both (4.4.14) and (4.4.15), as $\|\Pi\| \rightarrow 0$, the errors approach zero. However, before we let $\|\Pi\| \rightarrow 0$, we must first sum these equations over $j$, and the smaller we make $\|\Pi\|$, the more terms there are in the sum. When we sum both sides of (4.4.14), the errors accumulate, and although the error in each summand approaches zero as $\|\Pi\| \rightarrow 0$, the sum of the errors does not. When we use the more accurate approximation (4.4.15), this does not happen; the limit of the sum of the smaller errors is zero. We need the extra accuracy of (4.4.15) because the paths of Brownian motion are so volatile (i.e., they have nonzero quadratic variation). This extra term makes stochastic calculus different from ordinary calculus.

The Ito--Doeblin formula often simplifies the computation of Itô integrals. For example, with $f(x)=\frac{1}{2} x^{2}$, this formula says that

$$
\begin{aligned}
\frac{1}{2} W^{2}(T) & =f(W(T))-f(W(0)) \\
& =\int_{0}^{T} f^{\prime}(W(t)) d W(t)+\frac{1}{2} \int_{0}^{t} f^{\prime \prime}(W(t)) d t \\
& =\int_{0}^{T} W(t) d W(t)+\frac{1}{2} T
\end{aligned}
$$

Rearranging terms, we have formula (4.3.6) and have obtained it without going through the approximation of the integrand by simple processes as we did in Example 4.3.2.

\subsubsection{Formula for Itô Processes}

We extend the Itô-Doeblin formula to stochastic processes more general than Brownian motion. The processes for which we develop stochastic calculus are the Itô processes defined below. Almost all stochastic processes, except those that have jumps, are Itô processes.

Definition 4.4.3. Let $W(t), t \geq 0$, be a Brownian motion, and let $\mathcal{F}(t)$, $t \geq 0$, be an associated filtration. An Itô process is a stochastic process of the form

$$
X(t)=X(0)+\int_{0}^{t} \Delta(u) d W(u)+\int_{0}^{t} \Theta(u) d u
$$

where $X(0)$ is nonrandom and $\Delta(u)$ and $\Theta(u)$ are adapted stochastic processes. $^{2}$

In order to understand the volatility associated with Itô processes, we must determine the rate at which they accumulate quadratic variation.

Lemma 4.4.4. The quadratic variation of the Itô process (4.4.16) is

$$
[X, X](t)=\int_{0}^{t} \Delta^{2}(u) d u .
$$

Proof: We introduce the notation $I(t)=\int_{0}^{t} \Delta(u) d W(u), R(t)=\int_{0}^{t} \Theta(u) d u$. Both these processes are continuous in their upper limit of integration $t$. To

${ }^{2}$ We assume that $\mathbb{E} \int_{0}^{t} \Delta^{2}(u) d u$ and $\int_{0}^{t}|\Theta(u)| d u$ are finite for every $t>0$ so that the integrals on the right-hand side of (4.4.16) are defined and the Itô integral is a martingale. We shall always make such integrability assumptions, but we do not always explicitly state them. determine the quadratic variation of $X$ on $[0, t]$, we choose a partition $\Pi=$ $\left\{t_{0}, t_{1}, \ldots, t_{n}\right\}$ of $[0, t]$ (i.e., $0=t_{0}<t_{1}<\cdots<t_{n}=t$ ) and we write the sampled quadratic variation

$$
\begin{gathered}
\sum_{j=0}^{n-1}\left[X\left(t_{j+1}\right)-X\left(t_{j}\right)\right]^{2}=\sum_{j=0}^{n-1}\left[I\left(t_{j+1}\right)-I\left(t_{j}\right)\right]^{2}+\sum_{j=0}^{n-1}\left[R\left(t_{j+1}\right)-R\left(t_{j}\right)\right]^{2} \\
+2 \sum_{j=0}^{n-1}\left[I\left(t_{j+1}\right)-I\left(t_{j}\right)\right]\left[R\left(t_{j+1}\right)-R\left(t_{j}\right)\right] .
\end{gathered}
$$

As $\|\Pi\| \rightarrow 0$, the first term on the right-hand side, $\sum_{j=0}^{n-1}\left[I\left(t_{j+1}\right)-I\left(t_{j}\right)\right]^{2}$, converges to the quadratic variation of $I$ on $[0, t]$, which according to Theorem 4.3.1(vi) is $[I, I](t)=\int_{0}^{t} \Delta^{2}(u) d u$. The absolute value of the second term is bounded above by

$$
\begin{aligned}
& \max _{0 \leq k \leq n-1}\left|R\left(t_{k+1}\right)-R\left(t_{k}\right)\right| \cdot \sum_{j=0}^{n-1}\left|R\left(t_{j+1}\right)-R\left(t_{j}\right)\right| \\
&=\max _{0 \leq k \leq n-1}\left|R\left(t_{k+1}\right)-R\left(t_{k}\right)\right| \cdot \sum_{j=0}^{n-1}\left|\int_{t_{j}}^{t_{j+1}} \Theta(u) d u\right| \\
& \leq \max _{0 \leq k \leq n-1}\left|R\left(t_{k+1}\right)-R\left(t_{k}\right)\right| \cdot \sum_{j=0}^{n-1} \int_{t_{j}}^{t_{j+1}}|\Theta(u)| d u \\
&=\max _{0 \leq k \leq n-1}\left|R\left(t_{k+1}\right)-R\left(t_{k}\right)\right| \cdot \int_{0}^{t}|\Theta(u)| d u,
\end{aligned}
$$

and as $\|\Pi\| \rightarrow 0$, this has limit $0 \cdot \int_{0}^{t}|\Theta(u)| d u=0$ because $R(t)$ is continuous. The absolute value of the third term is bounded above by

$$
\begin{gathered}
2 \max _{0 \leq k \leq n-1}\left|I\left(t_{k+1}\right)-I\left(t_{k}\right)\right| \cdot \sum_{j=0}^{n-1}\left|R\left(t_{j+1}\right)-R\left(t_{j}\right)\right| \\
\quad \leq 2 \max _{0 \leq k \leq n-1}\left|I\left(t_{k+1}\right)-I\left(t_{k}\right)\right| \cdot \int_{0}^{t}|\Theta(u)| d u,
\end{gathered}
$$

and this has limit $0 \cdot \int_{0}^{t}|\Theta(u)|^{2} d u=0$ as $\|\Pi\| \rightarrow 0$ because $I(t)$ is continuous. We conclude that $[X, X](t)=[I, I](t)=\int_{0}^{t} \Delta^{2}(u) d u$.

The conclusion of Lemma 4.4.4 is most easily remembered by first writing (4.4.16) in the differential notation

$$
d X(t)=\Delta(t) d W(t)+\Theta(t) d t
$$

and then using the differential multiplication table (4.4.12) to compute 

$$
\begin{aligned}
d X(t) d X(t) & =\Delta^{2}(t) d W(t) d W(t)+2 \Delta(t) \Theta(t) d W(t) d t+\Theta^{2}(t) d t d t \\
& =\Delta^{2}(t) d t
\end{aligned}
$$

This says that, at each time $t$, the process $X$ is accumulating quadratic variation at rate $\Delta^{2}(t)$ per unit time, and hence the total quadratic variation accumulated on the time interval $[0, t]$ is $[X, X](t)=\int_{0}^{t} \Delta^{2}(u) d u$. This quadratic variation is solely due to the quadratic variation of the Itô integral $I(t)=\int_{0}^{t} \Delta(u) d W(u)$. The ordinary integral $R(t)=\int_{0}^{t} \Theta(u) d u$ has zero quadratic variation and thus contributes nothing to the quadratic variation of $X$.

Notice in this connection that having zero quadratic variation does not necessarily mean that $R(t)$ is nonrandom. Because $\Theta(u)$ can be random, $R(t)$ can also be random. However, $R(t)$ is not as volatile as $I(t)$. At each time $t$, we have a good estimate of the next increment of $R(t)$. For small time steps $h>0$

$$
R(t+h) \approx R(t)+\Theta(t) h,
$$

and we know both $R(t)$ and $\Theta(t)$ at time $t$. This is like investing in a money market account at a variable interest rate. At each time, we have a good estimate of the return over the near future because we know today's interest rate. Nonetheless, the return is random because the interest rate $(\theta$ in this analogy) can change. In contrast, $I$ is more volatile. At time $t$, one estimate of $I(t+h)$ is

$$
I(t+h) \approx I(t)+\Delta(t)(W(t+h)-W(t)),
$$

but we do not know $W(t+h)-W(t)$ at time $t$. In fact, $W(t+h)-W(t)$ is independent of the information available at time $t$. This is like investing in a stock.

So far we have discussed integrals with respect to time, such as $R(t)=$ $\int_{0}^{t} \Theta(u) d u$ appearing in (4.4.16) and Itô integrals (integrals with respect to Brownian motion) such as $I(t)=\int_{0}^{t} \Delta(u) d W(u)$, also appearing in (4.4.16). In addition, we shall need integrals with respect to Itô processes (i.e., integrals of the form $\int_{0}^{t} \Gamma(u) d X(u)$, where $\Gamma$ is some adapted process). We define such an integral by separating $d X(t)$ into a $d W(t)$ term and a $d t$ term as in (4.4.18).

Definition 4.4.5. Let $X(t), t \geq 0$, be an Itô process as described in Definition 4.4.3, and let $\Gamma(t), t \geq 0$, be an adapted process. We define the integral with respect to an Itô process ${ }^{3}$

$$
\int_{0}^{t} \Gamma(u) d X(u)=\int_{0}^{t} \Gamma(u) \Delta(u) d W(u)+\int_{0}^{t} \Gamma(u) \Theta(u) d u .
$$

We again work through the sketch of the proof of Theorem 4.4.1, but with the Itô process $X(t)$ replacing the Brownian motion $W(t)$. In place of (4.4.9), we now have

${ }^{3}$ We assume that $E \int_{0}^{t} \Gamma^{2}(u) \Delta^{2}(u) d u$ and $\int_{0}^{t}|\Gamma(u) \Theta(u)| d u$ are finite for every $t>0$ so that the integrals on the right-hand side of (4.4.20) are defined. 

$$
\begin{aligned}
& f(T, X(T))-f(0, X(0)) \\
& =\sum_{j=0}^{n-1} f_{t}\left(t_{j}, X\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)+\sum_{j=0}^{n-1} f_{x}\left(t_{j}, X\left(t_{j}\right)\right)\left(X\left(t_{j+1}\right)-X\left(t_{j}\right)\right) \\
& +\frac{1}{2} \sum_{j=0}^{n-1} f_{x x}\left(t_{j}, X\left(t_{j}\right)\right)\left(X\left(t_{j+1}\right)-X\left(t_{j}\right)\right)^{2} \\
& +\sum_{j=0}^{n-1} f_{t x}\left(t_{j}, X\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)\left(X\left(t_{j+1}\right)-X\left(t_{j}\right)\right) \\
& +\frac{1}{2} \sum_{j=0}^{n-1} f_{t t}\left(t_{j}, X\left(t_{j}\right)\right)\left(t_{j+1}-t_{j}\right)^{2}+\text { higher-order terms. }
\end{aligned}
$$

The last two sums on the right-hand side have zero limits as $\|\Pi\| \rightarrow 0$ for the same reasons the analogous terms have zero limits in the sketch of the proof of Theorem 4.4.1 (see (4.4.10) and (4.4.11)). The higher-order terms likewise have limit zero. The limit of the first term on the right-hand side of (4.4.21) is $\int_{0}^{T} f_{t}(t, X(t)) d t$. The limit of the second term is

$$
\int_{0}^{T} f_{x}(t, X(t)) d X(t)=\int_{0}^{T} f_{x}(t, X(t)) \Delta(t) d W(t)+\int_{0}^{T} f_{x}(t, X(t)) \Theta(t) d t
$$

Finally, the limit of the third term on the right-hand side of (4.4.19) is

$$
\frac{1}{2} \int_{0}^{T} f_{x x}(t, X(t)) d[X, X](t)=\frac{1}{2} \int_{0}^{T} f_{x x}(t, X(t)) \Delta^{2}(t) d t
$$

because the Itô process $X(t)$ accumulates quadratic variation at rate $\Delta^{2}(t)$ per unit time (Lemma 4.4.4). These considerations lead to the following generalization of Theorem 4.4.1.

Theorem 4.4.6 (Itô-Doeblin formula for an Itô process). Let $X(t)$, $t \geq 0$, be an Itô process as described in Definition 4.4.3, and let $f(t, x)$ be a function for which the partial derivatives $f_{t}(t, x), f_{x}(t, x)$, and $f_{x x}(t, x)$ are defined and continuous. Then, for every $T \geq 0$,

$$
\begin{aligned}
& f(T, X(T)) \\
& =f(0, X(0))+\int_{0}^{T} f_{t}(t, X(t)) d t+\int_{0}^{T} f_{x}(t, X(t)) d X(t) \\
& \quad+\frac{1}{2} \int_{0}^{T} f_{x x}(t, X(t)) d[X, X](t) \\
& =f(0, X(0))+\int_{0}^{T} f_{t}(t, X(t)) d t+\int_{0}^{T} f_{x}(t, X(t)) \Delta(t) d W(t) \\
& +\int_{0}^{T} f_{x}(t, X(t)) \Theta(t) d t+\frac{1}{2} \int_{0}^{T} f_{x x}(t, X(t)) \Delta^{2}(t) d t .
\end{aligned}
$$

Remark 4.4.7 (Summary of stochastic calculus). Theorem 4.4.6 is stated in mathematically precise language. Every term on the right-hand side has a solid definition, and in the end the right-hand side reduces to a sum of a nonrandom quantity $f(0, X(0))$, three ordinary (Lebesgue) integrals with respect to time, and an Itô integral.

However, it is easier to remember and use the result of this theorem if we recast it in differential notation. We may rewrite (4.4.22) as

$d f(t, X(t))=f_{t}(t, X(t)) d t+f_{x}(t, X(t)) d X(t)+\frac{1}{2} f_{x x}(t, X(t)) d X(t) d X(t)$.

The guiding principle here is that we write out the Taylor series expansion of $f(t, X(t))$ with respect to all its arguments, which in this case are $t$ and $X(t)$. We take this Taylor series expansion out to first order for every argument that has zero quadratic variation, which in this case is $t$, and we take the expansion out to second order for every argument that has nonzero quadratic variation, which in this case is $X(t)$.

We may reduce (4.4.23) to an expression that involves only $d t$ and $d W(t)$ by using the differential form (4.4.18) of the Itô process (i.e., $d X(t)=$ $\Delta(t) d W(t)+\Theta(t) d t)$ and the formula (4.4.19) for the rate at which $X(t)$ accumulates quadratic variation (i.e., $\left.d X(t) d X(t)=\Delta^{2}(t) d t\right)$. This is obtained by squaring the formula for $d X(t)$ and using the multiplication table (4.4.12). Making these substitutions in (4.4.23), we obtain

$$
\begin{aligned}
d f(t, X(t))=f_{t}(t, & X(t)) d t+f_{x}(t, X(t)) \Delta(t) d W(t) \\
& +f_{x}(t, X(t)) \Theta(t) d t+\frac{1}{2} f_{x x}(t, X(t)) \Delta^{2}(t) d t
\end{aligned}
$$

Itô calculus is little more than repeated use of this formula in a variety of situations.

\subsubsection{Examples}

We conclude this section with three examples illustrating Remark 4.4.7. Many more examples are developed in subsequent sections and in the exercises.

Example 4.4.8 (Generalized geometric Brownian motion). Let $W(t), t \geq 0$, be a Brownian motion, let $\mathcal{F}(t), t \geq 0$, be an associated filtration, and let $\alpha(t)$ and $\sigma(t)$ be adapted processes. Define the Itô process

$$
X(t)=\int_{0}^{t} \sigma(s) d W(s)+\int_{0}^{t}\left(\alpha(s)-\frac{1}{2} \sigma^{2}(s)\right) d s
$$

Then

$$
d X(t)=\sigma(t) d W(t)+\left(\alpha(t)-\frac{1}{2} \sigma^{2}(t)\right) d t
$$

and

$$
d X(t) d X(t)=\sigma^{2}(t) d W(t) d W(t)=\sigma^{2}(t) d t .
$$

Consider an asset price process given by

$$
S(t)=S(0) e^{X(t)}=S(0) \exp \left\{\int_{0}^{t} \sigma(s) d W(s)+\int_{0}^{t}\left(\alpha(s)-\frac{1}{2} \sigma^{2}(s)\right) d s\right\},
$$

where $S(0)$ is nonrandom and positive. We may write $S(t)=f(X(t))$, where $f(x)=S(0) e^{x}, f^{\prime}(x)=S(0) e^{x}$, and $f^{\prime \prime}(x)=S(0) e^{x}$. According to the ItôDoeblin formula

$$
\begin{aligned}
d S(t) & =d f(X(t)) \\
& =f^{\prime}(X(t)) d X(t)+\frac{1}{2} f^{\prime \prime}(X(t)) d X(t) d X(t) \\
& =S(0) e^{X(t)} d X(t)+\frac{1}{2} S(0) e^{X(t)} d X(t) d X(t) \\
& =S(t) d X(t)+\frac{1}{2} S(t) d X(t) d X(t) \\
& =\alpha(t) S(t) d t+\sigma(t) S(t) d W(t)
\end{aligned}
$$

The asset price $S(t)$ has instantaneous mean rate of return $\alpha(t)$ and volatility $\sigma(t)$. Both the instantaneous mean rate of return and the volatility are allowed to be time-varying and random.

This example includes all possible models of an asset price process that is always positive, has no jumps, and is driven by a single Brownian motion. Although the model is driven by a Brownian motion, the distribution of $S(t)$ does not need to be log-normal because $\alpha(t)$ and $\sigma(t)$ are allowed to be timevarying and random. If $\alpha$ and $\sigma$ are constant, we have the usual geometric Brownian motion model, and the distribution of $S(t)$ is log-normal.

In the case of constant $\alpha$ and $\sigma,(4.4 .26)$ becomes

$$
S(t)=S(0) \exp \left\{\sigma W(t)+\left(\alpha-\frac{1}{2} \sigma^{2}\right) t\right\} .
$$

One can incorrectly argue from this formula that since Brownian motion is a martingale (i.e., it has no overall tendency to rise or fall), the mean rate of return for $S(t)$ must be $\alpha-\frac{1}{2} \sigma^{2}$. The error in this argument is that although $W(t)$ is a martingale, $S(0) e^{\sigma W(t)}$ is not. The convexity of the function $e^{\sigma x}$ imparts an upward drift to $S(0) e^{\sigma W(t)}$. In order to correct for this, one must subtract $\frac{1}{2} \sigma^{2} t$ in the exponential; the process $S(0) \exp \left\{\sigma W(t)-\frac{1}{2} \sigma^{2} t\right\}$ is a martingale (see Theorem 3.6.1). If we now add $\alpha t$ in the exponential, we get $S(t)$, a process with mean rate of return $\alpha$.

The Itô-Doeblin formula automatically keeps track of these effects, even when $\alpha$ and $\sigma$ are time-varying and random. If $\alpha=0$, then (4.4.27) yields

$$
d S(t)=\sigma(t) S(t) d W(t) .
$$



\section{Integration of both sides yields}

$$
S(t)=S(0)+\int_{0}^{t} \sigma(s) S(s) d W(s) .
$$

The right-hand side is the nonrandom constant $S(0)$ plus an Itô integral, which is a martingale, and hence (in the case $\alpha=0$ )

$$
S(t)=S(0) \exp \left\{\int_{0}^{t} \sigma(s) d W(s)-\frac{1}{2} \int_{0}^{t} \sigma^{2}(s) d s\right\}
$$

is a martingale. In other words, the term $\sigma(t) S(t) d W(t)$ on the right-hand side of (4.4.27) contributes no drift, just pure volatility, to the asset price.

When $\alpha(t)$ is a nonzero random process, (4.4.27) shows that it plays the role of the mean rate of return. In the case of time-varying and random $\alpha(t)$, we will call this the instantaneous mean rate of return since it depends on the time (and the sample path) where it is evaluated.

The preceding example supplies the heart of the proof of the following theorem.

Theorem 4.4.9 (Itô integral of a deterministic integrand). Let $W(s)$, $s \geq 0$, be a Brownian motion, and let $\Delta(s)$ be a nonrandom function of time. Define $I(t)=\int_{0}^{t} \Delta(s) d W(s)$. For each $t \geq 0$, the random variable $I(t)$ is normally distributed with expected value zero and variance $\int_{0}^{t} \Delta^{2}(s) d s$.

Proof: The mean and variance of $I(t)$ are easy to determine. Since $I(t)$ is a martingale and $I(0)=0$, we must have $\mathbb{E} I(t)=I(0)=0$. Itô's isometry (Theorem 4.3.1(v)) implies that

$$
\operatorname{Var} I(t)=\mathbb{E} I^{2}(t)=\int_{0}^{t} \Delta^{2}(s) d s .
$$

We do not need to take the expected value of $\int_{0}^{t} \Delta^{2}(s) d s$ on the right-hand side of this formula because $\Delta(s)$ is not random.

The challenge is to show that $I(t)$ is normally distributed. We shall do this by establishing that $I(t)$ has the moment-generating function of a normal random variable with mean zero and variance $\int_{0}^{t} \Delta^{2}(s) d s$, which is (see $(3.2 .13))$

$$
\mathbb{E} e^{u I(t)}=\exp \left\{\frac{1}{2} u^{2} \int_{0}^{t} \Delta^{2}(s) d s\right\} \text { for all } u \in \mathbb{R} .
$$

Because $\Delta(s)$ is not random, (4.4.30) is equivalent to

$$
\mathbb{E} \exp \left\{u I(t)-\frac{1}{2} u^{2} \int_{0}^{t} \Delta^{2}(s) d s\right\}=1
$$

which may be rewritten as

$$
\mathbb{E} \exp \left\{\int_{0}^{t} u \Delta(s) d W(s)-\frac{1}{2} \int_{0}^{t}(u \Delta(s))^{2} d s\right\}=1 .
$$

But the process

$$
\exp \left\{\int_{0}^{t} u \Delta(s) d W(s)-\frac{1}{2} \int_{0}^{t}(u \Delta(s))^{2} d s\right\}
$$

is a martingale. Indeed, it is a generalized geometric Brownian motion with mean rate of return $\alpha=0$; see (4.4.29) with $\sigma(s)=u \Delta(s)$. Furthermore, this process takes the value 1 at $t=0$, and hence its expectation is always 1 . This gives us (4.4.31).

Note that (4.4.31) always holds, regardless of whether $\Delta(s)$ is random. However, we need to assume that $\Delta(s)$ is nonrandom in order to obtain the moment-generating function formula (4.4.30) from (4.4.31). When $\Delta(s)$ is random, there is no reason that the distribution of $\int_{0}^{t} \Delta(s) d W(s)$ should be normal.

Example 4.4.10 (Vasicek interest rate model). Let $W(t), t \geq 0$, be a Brownian motion. The Vasicek model for the interest rate process $R(t)$ is

$$
d R(t)=(\alpha-\beta R(t)) d t+\sigma d W(t),
$$

where $\alpha, \beta$, and $\sigma$ are positive constants. Equation (4.4.32) is an example of a stochastic differential equation. It defines a random process, $R(t)$ in this case, by giving a formula for its differential, and the formula involves the random process itself and the differential of a Brownian motion.

The solution to the stochastic differential equation (4.4.32) can be determined in closed form and is

$$
R(t)=e^{-\beta t} R(0)+\frac{\alpha}{\beta}\left(1-e^{-\beta t}\right)+\sigma e^{-\beta t} \int_{0}^{t} e^{\beta s} d W(s),
$$

a claim that we now verify. In particular, we compute the differential of the right-hand side of (4.4.33). To do this, we use the Itô-Doeblin formula with

$$
f(t, x)=e^{-\beta t} R(0)+\frac{\alpha}{\beta}\left(1-e^{-\beta t}\right)+\sigma e^{-\beta t} x
$$

and $X(t)=\int_{0}^{t} e^{\beta s} d W(s)$. Then the right-hand side of (4.4.33) is $f(t, X(t))$. The technique we are using is to separate the right-hand side into two parts: an ordinary function of two variables $t$ and $x$, which has no randomness in it, and an Itô process $X(t)$, which contains all the randomness. For the Itô-Doeblin formula, we shall need the following partial derivatives of $f(t, x)$ : 

$$
\begin{aligned}
f_{t}(t, x) & =-\beta e^{-\beta t} R(0)+\alpha e^{-\beta t}-\sigma \beta e^{-\beta t} x=\alpha-\beta f(t, x), \\
f_{x}(t, x) & =\sigma e^{-\beta t} \\
f_{x x}(t, x) & =0
\end{aligned}
$$

We shall also need the differential of $X(t)$, which is $d X(t)=e^{\beta t} d W(t)$. We shall not need $d X(t) d X(t)=e^{2 \beta t} d t$ because $f_{x x}(t, x)=0$. The Itô-Doeblin formula states that

$$
\begin{aligned}
& d f(t, X(t)) \\
& =f_{t}(t, X(t)) d t+f_{x}(t, X(t)) d X(t)+\frac{1}{2} f_{x x}(t, X(t)) d X(t) d X(t) \\
& =(\alpha-\beta f(t, X(t))) d t+\sigma d W(t)
\end{aligned}
$$

This shows that $f(t, X(t))$ satisfies the stochastic differential equation (4.4.32) that defines $R(t)$. Moreover, $f(0, X(0))=R(0)$. Because $f(t, X(t))$ satisfies the equation defining $R(t)$ and has the same initial condition as $R(t)$, it must be the case that $f(t, X(t))=R(t)$ for all $t \geq 0$.

Theorem 4.4.9 implies that the random variable $\int_{0}^{t} e^{\beta s} d W(s)$ appearing on the right-hand side of (4.4.33) is normally distributed with mean zero and variance

$$
\int_{0}^{t} e^{2 \beta s} d s=\frac{1}{2 \beta}\left(e^{2 \beta t}-1\right) .
$$

Therefore, $R(t)$ is normally distributed with mean $e^{-\beta t} R(0)+\frac{\alpha}{\beta}\left(1-e^{-\beta t}\right)$ and variance $\frac{\sigma^{2}}{2 \beta}\left(1-e^{-2 \beta t}\right)$. In particular, no matter how the parameters $\alpha>0$, $\beta>0$, and $\sigma>0$ are chosen, there is positive probability that $R(t)$ is negative, an undesirable property for an interest rate model.

The Vasicek model has the desirable property that the interest rate is mean-reverting. When $R(t)=\frac{\alpha}{\beta}$, the drift term (the $d t$ term) in (4.4.32) is zero. When $R(t)>\frac{\alpha}{\beta}$, this term is negative, which pushes $R(t)$ back toward $\frac{\alpha}{\beta}$. When $R(t)<\frac{\alpha}{\beta}$, this term is positive, which again pushes $R(t)$ back toward $\frac{\alpha}{\beta}$. If $R(0)=\frac{\alpha}{\beta}$, then $\mathbb{E} R(t)=\frac{\alpha}{\beta}$ for all $t \geq 0$. If $R(0) \neq \frac{\alpha}{\beta}$, then $\lim _{t \rightarrow \infty} \mathbb{E} R(t)=\frac{\alpha}{\beta}$.

Example 4.4.11 (Cox-Ingersoll-Ross (CIR) interest rate model). Let $W(t)$, $t \geq 0$, be a Brownian motion. The Cox-Ingersoll-Ross model for the interest rate process $R(t)$ is

$$
d R(t)=(\alpha-\beta R(t)) d t+\sigma \sqrt{R(t)} d W(t)
$$

where $\alpha, \beta$, and $\sigma$ are positive constants. Unlike the Vasicek equation (4.4.32), the CIR equation (4.4.34) does not have a closed-form solution. The advantage of (4.4.34) over the Vasicek model is that the interest rate in the CIR model does not become negative. If $R(t)$ reaches zero, the term multiplying $d W(t)$ vanishes and the positive drift term $\alpha d t$ in equation (4.4.34) drives the interest rate back into positive territory. Like the Vasicek model, the CIR model is mean-reverting.

Although one cannot derive a closed-form solution for (4.4.34), the distribution of $R(t)$ for each positive $t$ can be determined. That computation would take us too far afield. We instead content ourselves with the derivation of the expected value and variance of $R(t)$. To do this, we use the function $f(t, x)=e^{\beta t} x$ and the Itô-Doeblin formula to compute

$$
\begin{aligned}
& d\left(e^{\beta t} R(t)\right) \\
& =d f(t, R(t)) \\
& =f_{t}(t, R(t)) d t+f_{x}(t, R(t)) d R(t)+\frac{1}{2} f_{x x}(t, R(t)) d R(t) d R(t) \\
& =\beta e^{\beta t} R(t) d t+e^{\beta t}(\alpha-\beta R(t)) d t+e^{\beta t} \sigma \sqrt{R(t)} d W(t) \\
& =\alpha e^{\beta t} d t+\sigma e^{\beta t} \sqrt{R(t)} d W(t) .
\end{aligned}
$$

Integration of both sides of (4.4.35) yields

$$
\begin{aligned}
e^{\beta t} R(t) & =R(0)+\alpha \int_{0}^{t} e^{\beta u} d u+\sigma \int_{0}^{t} e^{\beta u} \sqrt{R(u)} d W(u) \\
& =R(0)+\frac{\alpha}{\beta}\left(e^{\beta t}-1\right)+\sigma \int_{0}^{t} e^{\beta u} \sqrt{R(u)} d W(u)
\end{aligned}
$$

Recalling that the expectation of an Itô integral is zero, we obtain

$$
e^{\beta t} \mathbb{E} R(t)=R(0)+\frac{\alpha}{\beta}\left(e^{\beta t}-1\right)
$$

or, equivalently,

$$
\mathbb{E} R(t)=e^{-\beta t} R(0)+\frac{\alpha}{\beta}\left(1-e^{-\beta t}\right) .
$$

This is the same expectation as in the Vasicek model.

To compute the variance of $R(t)$, we set $X(t)=e^{\beta t} R(t)$, for which we have already computed

$$
d X(t)=\alpha e^{\beta t} d t+\sigma e^{\beta t} \sqrt{R(t)} d W(t)=\alpha e^{\beta t} d t+\sigma e^{\frac{\beta t}{2}} \sqrt{X(t)} d W(t)
$$

and $\mathbb{E} X(t)=R(0)+\frac{\alpha}{\beta}\left(e^{\beta t}-1\right)$. According to the Itô-Doeblin formula (with $f(x)=x^{2}, f^{\prime}(x)=2 x$, and $\left.f^{\prime \prime}(x)=2\right)$

$$
\begin{aligned}
d\left(X^{2}(t)\right) & =2 X(t) d X(t)+d X(t) d X(t) \\
& =2 \alpha e^{\beta t} X(t) d t+2 \sigma e^{\frac{\beta t}{2}} X^{\frac{3}{2}}(t) d W(t)+\sigma^{2} e^{\beta t} X(t) d t .
\end{aligned}
$$

Integration of (4.4.37) yields 

$$
X^{2}(t)=X^{2}(0)+\left(2 \alpha+\sigma^{2}\right) \int_{0}^{t} e^{\beta u} X(u) d u+2 \sigma \int_{0}^{t} e^{\frac{\beta u}{2}} X^{\frac{3}{2}}(u) d W(u) .
$$

Taking expectations, using the fact that the expectation of an Itô integral is zero and the formula already derived for $\mathbb{E} X(t)$, we obtain

$$
\begin{aligned}
\mathbb{E} X^{2}(t) & =X^{2}(0)+\left(2 \alpha+\sigma^{2}\right) \int_{0}^{t} e^{\beta u} \mathbb{E} X(u) d u \\
& =R^{2}(0)+\left(2 \alpha+\sigma^{2}\right) \int_{0}^{t} e^{\beta u}\left(R(0)+\frac{\alpha}{\beta}\left(e^{\beta u}-1\right)\right) d u \\
& =R^{2}(0)+\frac{2 \alpha+\sigma^{2}}{\beta}\left(R(0)-\frac{\alpha}{\beta}\right)\left(e^{\beta t}-1\right)+\frac{2 \alpha+\sigma^{2}}{2 \beta} \cdot \frac{\alpha}{\beta}\left(e^{2 \beta t}-1\right) .
\end{aligned}
$$

Therefore,

$$
\begin{aligned}
\mathbb{E} R^{2}(t)= & e^{-2 \beta t} \mathbb{E} X^{2}(t) \\
= & e^{-2 \beta t} R^{2}(0)+\frac{2 \alpha+\sigma^{2}}{\beta}\left(R(0)-\frac{\alpha}{\beta}\right)\left(e^{-\beta t}-e^{-2 \beta t}\right) \\
& +\frac{\alpha\left(2 \alpha+\sigma^{2}\right)}{2 \beta^{2}}\left(1-e^{-2 \beta t}\right) .
\end{aligned}
$$

Finally,

$$
\begin{aligned}
\operatorname{Var}(R(t))= & \mathbb{E} R^{2}(t)-(\mathbb{E} R(t))^{2} \\
= & e^{-2 \beta t} R^{2}(0)+\frac{2 \alpha+\sigma^{2}}{\beta}\left(R(0)-\frac{\alpha}{\beta}\right)\left(e^{-\beta t}-e^{-2 \beta t}\right) \\
& +\frac{\alpha\left(2 \alpha+\sigma^{2}\right)}{2 \beta^{2}}\left(1-e^{-2 \beta t}\right)-e^{-2 \beta t} R^{2}(0) \\
& \quad-\frac{2 \alpha}{\beta} R(0)\left(e^{-\beta t}-e^{-2 \beta t}\right)-\frac{\alpha^{2}}{\beta^{2}}\left(1-e^{-\beta t}\right)^{2} \\
= & \frac{\sigma^{2}}{\beta} R(0)\left(e^{-\beta t}-e^{-2 \beta t}\right)+\frac{\alpha \sigma^{2}}{2 \beta^{2}}\left(1-2 e^{-\beta t}+e^{-2 \beta t}\right) .
\end{aligned}
$$

In particular,

$$
\lim _{t \rightarrow \infty} \operatorname{Var}(R(t))=\frac{\alpha \sigma^{2}}{2 \beta^{2}}
$$

\subsection{Black-Scholes-Merton Equation}

The addition of Merton's name to what has traditionally been called the Black-Scholes equation is explained in the Notes, Section $4.9$. In this section, we derive the Black-Scholes-Merton partial differential equation for the price of an option on an asset modeled as a geometric Brownian motion. The idea behind this derivation is the same as in the binomial model of Chapter 1 of Volume I, which is to determine the initial capital required to perfectly hedge a short position in the option.

\subsubsection{Evolution of Portfolio Value}

Consider an agent who at each time $t$ has a portfolio valued at $X(t)$. This portfolio invests in a money market account paying a constant rate of interest $r$ and in a stock modeled by the geometric Brownian motion

$$
d S(t)=\alpha S(t) d t+\sigma S(t) d W(t)
$$

Suppose at each time $t$, the investor holds $\Delta(t)$ shares of stock. The position $\Delta(t)$ can be random but must be adapted to the filtration associated with the Brownian motion $W(t), t \geq 0$. The remainder of the portfolio value, $X(t)-\Delta(t) S(t)$, is invested in the money market account.

The differential $d X(t)$ for the investor's portfolio value at each time $t$ is due to two factors, the capital gain $\Delta(t) d S(t)$ on the stock position and the interest earnings $r(X(t)-\Delta(t) S(t)) d t$ on the cash position. In other words,

$$
\begin{aligned}
d X(t) & =\Delta(t) d S(t)+r(X(t)-\Delta(t) S(t)) d t \\
& =\Delta(t)(\alpha S(t) d t+\sigma S(t) d W(t))+r(X(t)-\Delta(t) S(t)) d t \\
& =r X(t) d t+\Delta(t)(\alpha-r) S(t) d t+\Delta(t) \sigma S(t) d W(t) .
\end{aligned}
$$

The three terms appearing in the last line of (4.5.2) can be understood as follows:

(i) an average underlying rate of return $r$ on the portfolio, which is reflected by the term $r X(t) d t$,

(ii) a risk premium $\alpha-r$ for investing in the stock, which is reflected by the term $\Delta(t)(\alpha-r) S(t) d t$, and

(iii) a volatility term proportional to the size of the stock investment, which is the term $\Delta(t) \sigma S(t) d W(t)$.

The discrete-time analogue of equation (4.5.2) appears in Chapter 1 of Volume I as (1.2.12):

$$
X_{n+1}=\Delta_{n} S_{n+1}+(1+r)\left(X_{n}-\Delta_{n} S_{n}\right) .
$$

We may rearrange terms in this equation to obtain

$$
X_{n+1}-X_{n}=\Delta_{n}\left(S_{n+1}-S_{n}\right)+r\left(X_{n}-\Delta_{n} S_{n}\right),
$$

which is analogous to the first line of (4.5.2), except in (4.5.3) time steps forward one unit at a time, whereas in (4.5.2) time moves forward continuously. See Exercise $4.10$ for additional discussion of the rationale for equation (4.5.2) in option pricing.

We shall of ten consider the discounted stock price $e^{-r t} S(t)$ and the discounted portfolio value of an agent, $e^{-r t} X(t)$. According to the Itô-Doeblin formula with $f(t, x)=e^{-r t} x$, the differential of the discounted stock price is

$$
\begin{aligned}
& d\left(e^{-r t} S(t)\right) \\
& =d f(t, S(t)) \\
& =f_{t}(t, S(t)) d t+f_{x}(t, S(t)) d S(t)+\frac{1}{2} f_{x x}(t, S(t)) d S(t) d S(t) \\
& =-r e^{-r t} S(t) d t+e^{-r t} d S(t) \\
& =(\alpha-r) e^{-r t} S(t) d t+\sigma e^{-r t} S(t) d W(t)
\end{aligned}
$$

and the differential of the discounted portfolio value is

$$
\begin{aligned}
& d\left(e^{-r t} X(t)\right) \\
& =d f(t, X(t)) \\
& =f_{t}(t, X(t)) d t+f_{x}(t, X(t)) d X(t)+\frac{1}{2} f_{x x}(t, X(t)) d X(t) d X(t) \\
& =-r e^{-r t} X(t) d t+e^{-r t} d X(t) \\
& =\Delta(t)(\alpha-r) e^{-r t} S(t) d t+\Delta(t) \sigma e^{-r t} S(t) d W(t) \\
& =\Delta(t) d\left(e^{-r t} S(t)\right) .
\end{aligned}
$$

Discounting the stock price reduces the mean rate of return from $\alpha$, the term multiplying $S(t) d t$ in (4.5.1), to $\alpha-r$, the term multiplying $e^{-r t} S(t) d t$ in (4.5.4). Discounting the portfolio value removes the underlying rate of return $r$; compare the last line of (4.5.2) to the next-to-last line of (4.5.5). The last line of (4.5.5) shows that change in the discounted portfolio value is solely due to change in the discounted stock price.

\subsubsection{Evolution of Option Value}

Consider a European call option that pays $(S(T)-K)^{+}$at time $T$. The strike price $K$ is some nonnegative constant. Black, Scholes, and Merton argued that the value of this call at any time should depend on the time (more precisely, on the time to expiration) and on the value of the stock price at that time, and of course it should also depend on the model parameters $r$ and $\sigma$ and the contractual strike price $K$. Only two of these quantities, time and stock price, are variable. Following this reasoning, we let $c(t, x)$ denote the value of the call at time $t$ if the stock price at that time is $S(t)=x$. There is nothing random about the function $c(t, x)$. However, the value of the option is random; it is the stochastic process $c(t, S(t))$ obtained by replacing the dummy variable $x$ by the random stock price $S(t)$ in this function. At the initial time, we do not know the future stock prices $S(t)$ and hence do not know the future option values $c(t, S(t))$. Our goal is to determine the function $c(t, x)$ so we at least have a formula for the future option values in terms of the future stock prices.

We begin by computing the differential of $c(t, S(t))$. According to the ItôDoeblin formula, it is

$$
\begin{aligned}
& d c(t, S(t)) \\
& =c_{t}(t, S(t)) d t+c_{x}(t, S(t)) d S(t)+\frac{1}{2} c_{x x}(t, S(t)) d S(t) d S(t) \\
& =c_{t}(t, S(t)) d t+c_{x}(t, S(t))(\alpha S(t) d t+\sigma S(t) d W(t)) \\
& \quad \quad+\frac{1}{2} c_{x x}(t, S(t)) \sigma^{2} S^{2}(t) d t \\
& =\left[c_{t}(t, S(t))+\alpha S(t) c_{x}(t, S(t))+\frac{1}{2} \sigma^{2} S^{2}(t) c_{x x}(t, S(t))\right] d t \\
& \quad+\sigma S(t) c_{x}(t, S(t)) d W(t) .
\end{aligned}
$$

We next compute the differential of the discounted option price $e^{-r t} c(t, S(t))$. Let $f(t, x)=e^{-r t} x$. According to the Itô-Doeblin formula,

$$
\begin{aligned}
d\left(e^{-r t}\right. & c(t, S(t))) \\
= & d f(t, c(t, S(t))) \\
= & f_{t}(t, c(t, S(t))) d t+f_{x}(t, c(t, S(t))) d c(t, S(t)) \\
& \quad+\frac{1}{2} f_{x x}(t, c(t, S(t))) d c(t, S(t)) d c(t, S(t)) \\
= & -r e^{-r t} c(t, S(t)) d t+e^{-r t} d c(t, S(t)) \\
= & e^{-r t}\left[-r c(t, S(t))+c_{t}(t, S(t))+\alpha S(t) c_{x}(t, S(t))\right. \\
& \left.\quad+\frac{1}{2} \sigma^{2} S^{2}(t) c_{x x}(t, S(t))\right] d t+e^{-r t} \sigma S(t) c_{x}(t, S(t)) d W(t) .
\end{aligned}
$$

\subsubsection{Equating the Evolutions}

A (short option) hedging portfolio starts with some initial capital $X(0)$ and invests in the stock and money market account so that the portfolio value $X(t)$ at each time $t \in[0, T]$ agrees with $c(t, S(t))$. This happens if and only if $e^{-r t} X(t)=e^{-r t} c(t, S(t))$ for all $t$. One way to ensure this equality is to make sure that

$$
d\left(e^{-r t} X(t)\right)=d\left(e^{-r t} c(t, S(t))\right) \text { for all } t \in[0, T)
$$

and $X(0)=c(0, S(0))$. Integration of (4.5.8) from 0 to $t$ then yields

$$
e^{-r t} X(t)-X(0)=e^{-r t} c(t, S(t))-c(0, S(0)) \text { for all } t \in[0, T) .
$$

If $X(0)=c(0, S(0))$, then we can cancel this term in (4.5.9) and get the desired equality. Comparing (4.5.5) and (4.5.7), we see that (4.5.8) holds if and only if

$$
\begin{aligned}
& \Delta(t)(\alpha-r) S(t) d t+\Delta(t) \sigma S(t) d W(t) \\
& =\left[-r c(t, S(t))+c_{t}(t, S(t))+\alpha S(t) c_{x}(t, S(t))+\frac{1}{2} \sigma^{2} S^{2}(t) c_{x x}(t, S(t))\right] d t \\
& \quad+\sigma S(t) c_{x}(t, S(t)) d W(t)
\end{aligned}
$$

We examine what is required in order for (4.5.10) to hold.

We first equate the $d W(t)$ terms in (4.5.10), which gives

$$
\Delta(t)=c_{x}(t, S(t)) \text { for all } t \in[0, T) .
$$

This is called the delta-hedging rule. At each time $t$ prior to expiration, the number of shares held by the hedge of the short option position is the partial derivative with respect to the stock price of the option value at that time. This quantity, $c_{x}(t, S(t))$, is called the delta of the option.

We next equate the $d t$ terms in (4.5.10), using (4.5.11), to obtain

$$
\begin{aligned}
& (\alpha-r) S(t) c_{x}(t, S(t)) \\
& =-r c(t, S(t))+c_{t}(t, S(t))+\alpha S(t) c_{x}(t, S(t))+\frac{1}{2} \sigma^{2} S^{2}(t) c_{x x}(t, S(t)) \\
& \text { for all } t \in[0, T)
\end{aligned}
$$

The term $\alpha S(t) c_{x}(t, S(t))$ appears on both sides of (4.5.12), and after canceling it, we obtain

$$
\begin{aligned}
r c(t, S(t))=c_{t}(t, S(t))+r S(t) c_{x}(t, S(t))+ & \frac{1}{2} \sigma^{2} S^{2}(t) c_{x x}(t, S(t)) \\
& \text { for all } t \in[0, T)
\end{aligned}
$$

In conclusion, we should seek a continuous function $c(t, x)$ that is a solution to the Black-Scholes-Merton partial differential equation

$$
c_{t}(t, x)+r x c_{x}(t, x)+\frac{1}{2} \sigma^{2} x^{2} c_{x x}(t, x)=r c(t, x) \text { for all } t \in[0, T), x \geq 0
$$

and that satisfies the terminal condition

$$
c(T, x)=(x-K)^{+} .
$$

Suppose we have found this function. If an investor starts with initial capital $X(0)=c(0, S(0))$ and uses the hedge $\Delta(t)=c_{x}(t, S(t))$, then (4.5.10) will hold for all $t \in[0, T)$. Indeed, the $d W(t)$ terms on the left and right sides of (4.5.10) agree because $\Delta(t)=c_{x}(t, S(t))$, and the $d t$ terms agree because (4.5.14) guarantees (4.5.13). Equality in (4.5.10) gives us (4.5.9). Canceling $X(0)=c(0, S(0))$ and $e^{-r t}$ in this equation, we see that $X(t)=c(t, S(t))$ for all $t \in[0, T)$. Taking the limit as $t \uparrow T$ and using the fact that both $X(t)$ and $c(t, S(t))$ are continuous, we conclude that $X(T)=c(T, S(T))=(S(T)-K)^{+}$. This means that the short position has been successfully hedged. No matter which of its possible paths the stock price follows, when the option expires, the agent hedging the short position has a portfolio whose value agrees with the option payoff.

\subsubsection{Solution to the Black-Scholes-Merton Equation}

The Black-Scholes-Merton equation (4.5.14) does not involve probability. It is a partial differential equation, and the arguments $t$ and $x$ are dummy variables, not random variables. One can solve it by partial differential equation methods. In this section, however, rather than showing how to solve the equation, we shall simply present the solution and check that it works. In Subsection 5.2.5, we present a derivation of this solution based on probability theory.

We want the Black-Scholes-Merton equation to hold for all $x \geq 0$ and $t \in[0, T)$ so that (4.5.14) will hold regardless of which of its possible paths the stock price follows. If the initial stock price is positive, then the stock price is always positive, and it can take any positive value. If the initial stock price is zero, then subsequent stock prices are all zero. We cover both of these cases by asking (4.5.14) to hold for all $x \geq 0$. We do not need (4.5.14) to hold at $t=T$, although we need the function $c(t, x)$ to be continuous at $t=T$. If the hedge works at all times strictly prior to $T$, it also works at time $T$ because of continuity.

Equation (4.5.14) is a partial differential equation of the type called backward parabolic. For such an equation, in addition to the terminal condition (4.5.15), one needs boundary conditions at $x=0$ and $x=\infty$ in order to determine the solution. The boundary condition at $x=0$ is obtained by substituting $x=0$ into (4.5.14), which then becomes

$$
c_{t}(t, 0)=r c(t, 0) .
$$

This is an ordinary differential equation for the function $c(t, 0)$ of $t$, and the solution is

$$
c(t, 0)=e^{r t} c(0,0) .
$$

Substituting $t=T$ into this equation and using the fact that $c(T, 0)=(0-$ $K)^{+}=0$, we see that $c(0,0)=0$ and hence

$$
c(t, 0)=0 \text { for all } t \in[0, T] .
$$

This is the boundary condition at $x=0$.

As $x \rightarrow \infty$, the function $c(t, x)$ grows without bound. In such a case, we give the boundary condition at $x=\infty$ by specifying the rate of growth. One way to specify a boundary condition at $x=\infty$ for the European call is

$$
\lim _{x \rightarrow \infty}\left[c(t, x)-\left(x-e^{-r(T-t)} K\right)\right]=0 \text { for all } t \in[0, T] .
$$

In particular, $c(t, x)$ grows at the same rate as $x$ as $x \rightarrow \infty$. Recall that $c(t, x)$ is the value at time $t$ of a call on a stock whose price at time $t$ is $x$. For large $x$, this call is deep in the money and very likely to end in the money. In this case, the price of the call is almost as much as the price of the forward contract discussed in Subsection 4.5.6 below (see (4.5.26)). This is the assertion of (4.5.18).

The solution to the Black-Scholes-Merton equation (4.5.14) with terminal condition (4.5.15) and boundary conditions (4.5.17) and (4.5.18) is

$$
c(t, x)=x N\left(d_{+}(T-t, x)\right)-K e^{-r(T-t)} N\left(d_{-}(T-t, x)\right), 0 \leq t<T, x>0,
$$

where

$$
d_{\pm}(\tau, x)=\frac{1}{\sigma \sqrt{\tau}}\left[\log \frac{x}{K}+\left(r \pm \frac{\sigma^{2}}{2}\right) \tau\right]
$$

and $N$ is the cumulative standard normal distribution

$$
N(y)=\frac{1}{\sqrt{2 \pi}} \int_{-\infty}^{y} e^{-\frac{z^{2}}{2}} d z=\frac{1}{\sqrt{2 \pi}} \int_{-y}^{\infty} e^{-\frac{z^{2}}{2}} d z .
$$

We shall sometimes use the notation

$$
\operatorname{BSM}(\tau, x ; K, r, \sigma)=x N\left(d_{+}(\tau, x)\right)-K e^{-r \tau} N\left(d_{-}(\tau, x)\right),
$$

and call $\operatorname{BSM}(\tau, x ; K, r, \sigma)$ the Black-Scholes-Merton function. In this formula, $\tau$ and $x$ denote the time to expiration and the current stock price, respectively. The parameters $K, r$, and $\sigma$ are the strike price, the interest rate, and the stock volatility, respectively.

Formula (4.5.19) does not define $c(t, x)$ when $t=T$ (because then $\tau=T-$ $t=0$ and this appears in the denominator in (4.5.20)), nor does it define $c(t, x)$ when $x=0$ (because $\log x$ appears in (4.5.20), and $\log 0$ is not a real number). However, (4.5.19) defines $c(t, x)$ in such a way that $\lim _{t \rightarrow T} c(t, x)=(x-K)^{+}$ and $\lim _{x \downarrow 0} c(t, x)=0$. Verification of all of these claims is given as Exercise $4.9$.

\subsubsection{The Greeks}

The derivatives of the function $c(t, x)$ of (4.5.19) with respect to various variables are called the Greeks. Two of these are derived in Exercise 4.9, namely delta, which is

$$
c_{x}(t, x)=N\left(d_{+}(T-t, x)\right),
$$

and theta, which is

$$
c_{t}(t, x)=-r K e^{-r(T-t)} N\left(d_{-}(T-t, x)\right)-\frac{\sigma x}{2 \sqrt{T-t}} N^{\prime}\left(d_{+}(T-t, x)\right) .
$$

Because both $N$ and $N^{\prime}$ are always positive, delta is always positive and theta is always negative. Another of the Greeks is gamma, which is 

$$
c_{x x}(t, x)=N^{\prime}\left(d_{+}(T-t, x)\right) \frac{\partial}{\partial x} d_{+}(T-t, x)=\frac{1}{\sigma x \sqrt{T-t}} N^{\prime}\left(d_{+}(T-t, x)\right)
$$

Like delta, gamma is always positive.

In order to simplify notation in the following discussion, we sometimes suppress the arguments $(t, x)$ of $c(t, x)$ and $(T-t, x)$ of $d_{\pm}(T-t, x)$. If at time $t$ the stock price is $x$, then the short option hedge of (4.5.11) calls for holding $c_{x}(t, x)$ shares of stock, a position whose value is $x c_{x}=x N\left(d_{+}\right)$. The hedging portfolio value is $c=x N\left(d_{+}\right)-K e^{-r(T-t)} N\left(d_{-}\right)$, and since $x c_{x}(t, x)$ of this value is invested in stock, the amount invested in the money market must be

$$
c(t, x)-x c_{x}(t, x)=-K e^{-r(T-t)} N\left(d_{-}\right),
$$

a negative number. To hedge a short position in a call option, one must borrow money. To hedge a long position in a call option, one does the opposite. In other words, to hedge a long call position one should hold $-c_{x}$ shares of stock (i.e., have a short position in stock) and invest $K e^{-r(T-t)} N\left(d_{-}\right)$in the money market account.

Because delta and gamma are positive, for fixed $t$, the function $c(t, x)$ is increasing and convex in the variable $x$, as shown in Figure 4.5.1. Suppose at time $t$ the stock price is $x_{1}$ and we wish to take a long position in the option and hedge it. We do this by purchasing the option for $c\left(t, x_{1}\right)$, shorting $c_{x}\left(t, x_{1}\right)$ shares of stock, which generates income $x_{1} c_{x}\left(t, x_{1}\right)$, and investing the difference,

$$
M=x_{1} c_{x}\left(t, x_{1}\right)-c\left(t, x_{1}\right),
$$

in the money market account. We wish to consider the sensitivity to stock price changes of the portfolio that has these three components: long option, short stock, and long money market account. The initial portfolio value

$$
c\left(t, x_{1}\right)-x_{1} c_{x}\left(t, x_{1}\right)+M
$$

is zero at the moment $t$ when we set up these positions.

If the stock price were to instantaneously fall to $x_{0}$ as shown in Figure 4.5.1 and we do not change our positions in the stock or money market account, then the value of the option we hold would fall to $c\left(t, x_{0}\right)$ and the liability due to our short position in stock would decrease to $x_{0} c_{x}\left(t, x_{1}\right)$. Our total portfolio value, including $M$ in the money market account, would be

$$
c\left(t, x_{0}\right)-x_{0} c_{x}\left(t, x_{1}\right)+M=c\left(t, x_{0}\right)-c_{x}\left(t, x_{1}\right)\left(x_{0}-x_{1}\right)-c\left(t, x_{1}\right) .
$$

This is the difference at $x_{0}$ between the curve $y=c(t, x)$ and the straight line $y=c_{x}\left(t, x_{1}\right)\left(x-x_{1}\right)+c\left(t, x_{1}\right)$ in Figure 4.5.1. Because this difference is positive, our portfolio benefits from an instantaneous drop in the stock price.

On the other hand, if the stock price were to instantaneously rise to $x_{2}$ and we do not change our positions in the stock or money market account, then the value of the option would rise to $c\left(t, x_{2}\right)$ and the liability due to our short position in stock would increase to $x_{2} c_{x}\left(t, x_{1}\right)$. Our total portfolio value, including $M$ in the money market account, would be

$$
c\left(t, x_{2}\right)-x_{2} c_{x}\left(t, x_{1}\right)+M=c\left(t, x_{2}\right)-c_{x}\left(t, x_{1}\right)\left(x_{2}-x_{1}\right)-c\left(t, x_{1}\right) .
$$

This is the difference at $x_{2}$ between the curve $y=c(t, x)$ and the straight line $y=c_{x}\left(t, x_{1}\right)\left(x-x_{1}\right)+c\left(t, x_{1}\right)$ in Figure 4.5.1. This difference is positive, so our portfolio benefits from an instantaneous rise in the stock price.

![](https://cdn.mathpix.com/cropped/2023_02_09_f310c50a6371948e10fbg-181.jpg?height=479&width=1148&top_left_y=486&top_left_x=81)

Fig. 4.5.1. Delta-neutral position.

The portfolio we have set up is said to be delta-neutral and long gamma. The portfolio is long gamma because it benefits from the convexity of $c(t, x)$ as described above. If there is an instantaneous rise or an instantaneous fall in the stock price, the value of the portfolio increases. A long gamma portfolio is profitable in times of high stock volatility.

"Delta-neutral" refers to the fact that the line in Figure 4.5.1 is tangent to the curve $y=c(t, x)$. Therefore, when the stock price makes a small move, the change of portfolio value due to the corresponding change in option price is nearly offset by the change in the value of our short position in the stock. The straight line is a good approximation to the option price for small stock price moves. If the straight line were steeper than the option price curve at the starting point $x_{1}$, then we would be short delta; an upward move in the stock price would hurt the portfolio because the liability from the short position in stock would rise faster than the value of the option. On the other hand, a downward move would increase the portfolio value because the option price would fall more slowly than the rate of decrease in the liability from the short stock position. Unless a trader has a view on the market, he tries to set up portfolios that are delta-neutral. If he expects high volatility, he would at the same time try to choose the portfolio to be long gamma.

The portfolio described above may at first appear to offer an arbitrage opportunity. When we let time move forward, not only does the long gamma position offer an opportunity for profit, but the positive investment in the money market account enhances this opportunity. The drawback is that theta, the derivative of $c(t, x)$ with respect to time, is negative. As we move forward in time, the curve $y=c(t, x)$ is shifting downward. Figure 4.5.1 is misleading because it is drawn with $t$ fixed. In principle, the portfolio can lose money because the curve $c(t, x)$ shifts downward more rapidly than the money market investment and the long gamma position generate income. The essence of the hedging argument in Subsection 4.5.3 is that if the stock really is a geometric Brownian motion and we have determined the right value of the volatility $\sigma$, then so long as we continuously rebalance our portfolio, all these effects exactly cancel!

Of course, assets are not really geometric Brownian motions with constant volatility, but the argument above gives a good first approximation to reality. It also highlights volatility as the key parameter. In fact, the mean rate of return $\alpha$ of the stock does not appear in the Black-Scholes-Merton equation (4.5.14). From the point of view of no-arbitrage pricing, it is irrelevant how likely the stock is to go up or down because a delta-neutral position is a hedge against both possibilities. What matters is how much volatility the stock has, for we need to know the amount of profit that can be made from the long gamma position. The more volatile stocks offer more opportunity for profit from the portfolio that hedges a long call position with a short stock position, and hence the call is more expensive. The derivative of the option price with respect to the volatility $\sigma$ is called vega, and it is positive. As volatility increases, so do option prices in the Black-Scholes-Merton model.

\subsubsection{Put-Call Parity}

A forward contract with delivery price $K$ obligates its holder to buy one share of the stock at expiration time $T$ in exchange for payment $K$. At expiration, the value of the forward contract is $S(T)-K$. Let $f(t, x)$ denote the value of the forward contract at earlier times $t \in[0, T]$ if the stock price at time $t$ is $S(t)=x$.

We argue that the value of a forward contract is given by

$$
f(t, x)=x-e^{-r(T-t)} K .
$$

If an agent sells this forward contract at time zero for $f(t, S(0))=S(0)-$ $e^{-r T} K$, he can set up a static hedge, a hedge that does not trade except at the initial time, in order to protect himself. Specifically, the agent should purchase one share of stock. Since he has initial capital $S(0)-e^{-r T} K$ from the sale of the forward contract, this requires that he borrow $e^{-r T} K$ from the money market account. The agent makes no further trades. At expiration of the forward contract, he owns one share of stock and his debt to the money market account has grown to $K$, so his portfolio value is $S(T)-K$, exactly the value of the forward contract. Because the agent has been able to replicate the payoff of the forward contract with a portfolio whose value at each time $t$ is $S(t)-e^{-r(T-t)} K$, this must be the value at each time of the forward contract. This is $f(t, S(t))$, where $f(t, x)$ is defined by (4.5.26).

The forward price of a stock at time $t$ is defined to be the value of $K$ that causes the forward contract at time $t$ to have value zero (i.e., it is the value of $K$ that satisfies the equation $\left.S(t)-e^{-r(T-t)} K=0\right)$. Hence, we see that in a model with a constant interest rate, the forward price at time $t$ is

$$
\text { For }(t)=e^{r(T-t)} S(t) .
$$

Note that the forward price is not the price (or value) of a forward contract. For $0 \leq t \leq T$, the forward price at time $t$ is the price one can lock in at time $t$ for the purchase of one share of stock at time $T$, paying the price (settling) at time $T$. No money changes hands at the time the price is locked in.

Let us consider this situation at time $t=0$. At that time, one can lock in a price For $(0)=e^{r T} S(0)$ for purchase of the stock at time $T$. Let us do this, which means we set $K=e^{r T} S(0)$ in (4.5.26). The value of this forward contract is zero at time $t=0$, but as soon as time begins to move forward, the value of the forward contract changes. Indeed, its value at time $t$ is

$$
f(t, S(t))=S(t)-e^{r t} S(0) .
$$

Finally, let us consider a European put, which pays off $(K-S(T))^{+}$at time $T$. We observe that for any number $x$, the equation

$$
x-K=(x-K)^{+}-(K-x)^{+}
$$

holds. Indeed, if $x \geq K$, then $(x-K)^{+}=x-K$ and $(K-x)^{+}=0$. On the other hand, if $x \leq K$, then $(x-K)^{+}=0$ and $-(K-x)^{+}=-(K-x)=x-K$. In either case, the right-hand side of (4.5.28) equals the left-hand side. We denote by $p(t, x)$ the value of the European put at time $t$ if the time- $t$ stock price is $S(t)=x$. Similarly, we denote by $c(t, x)$ the value of the European call expiring at time $T$ with strike price $K$ and by $f(t, x)$ the value of the forward contract for the purchase of one share of stock at time $T$ in exchange for payment $K$. Equation (4.5.28) implies

$$
f(T, S(T))=c(T, S(T))-p(T, S(T)) ;
$$

the payoff of the forward contract agrees with the payoff of a portfolio that is long one call and short one put. Since the value at time $T$ of the forward contract agrees with the value of the portfolio that is long one call and short one put, these values must agree at all previous times:

$$
f(t, x)=c(t, x)-p(t, x), \quad x \geq 0,0 \leq t \leq T .
$$

If this were not the case, one could at some time $t$ either sell or buy the portfolio that is long the forward, short the call, and long the put, realizing an instant profit, and have no liability upon expiration of the contracts. The relationship (4.5.29) is called put-call parity. Note that we have derived the put-call parity formula (4.5.29) without appealing to the Black-Scholes-Merton model of a geometric Brownian motion for the stock price and a constant interest rate. Indeed, without any assumptions on the prices except sufficient liquidity that permits one to form the portfolio that is long one call and short one put, we have put-call parity. If we make the assumption of a constant interest rate $r$, then $f(t, x)$ is given by (4.5.26). If we make the additional assumption that the stock is a geometric Brownian motion with constant volatility $\sigma>0$, then we have also the BlackScholes-Merton call formula (4.5.19). We can then solve (4.5.29) to obtain the Black-Scholes-Merton put formula

$$
\begin{aligned}
p(t, x) & =x\left(N\left(d_{+}(T-t, x)\right)-1\right)-K e^{-r(T-t)}\left(N\left(d_{-}(T-t, x)\right)-1\right) \\
& =K e^{-r(T-t)} N\left(-d_{-}(T-t, x)\right)-x N\left(-d_{+}(T-t, x)\right),
\end{aligned}
$$

where $d_{\pm}(T-t, x)$ is given by (4.5.20).

\subsection{Multivariable Stochastic Calculus}

\subsubsection{Multiple Brownian Motions}

Definition 4.6.1. $A d$-dimensional Brownian motion is a process

$$
W(t)=\left(W_{1}(t), \ldots, W_{d}(t)\right)
$$

with the following properties.

(i) Each $W_{i}(t)$ is a one-dimensional Brownian motion.

(ii) If $i \neq j$, then the processes $W_{i}(t)$ and $W_{j}(t)$ are independent.

Associated with a d-dimensional Brownian motion, we have a filtration $\mathcal{F}(t)$, $t \geq 0$, such that the following holds.

(üi) (Information accumulates) For $0 \leq s<t$, every set in $\mathcal{F}(s)$ is also in $\mathcal{F}(t)$

(iv) (Adaptivity) For each $t \geq 0$, the random vector $W(t)$ is $\mathcal{F}(t)$-measurable.

(v) (Independence of future increments) For $0 \leq t<u$, the vector of increments $W(u)-W(t)$ is independent of $\mathcal{F}(t)$.

Although we have defined a multidimensional Brownian motion to be a vector of independent one-dimensional Brownian motions, we shall see in Example 4.6.6 how to build correlated Brownian motions from this.

Because each component $W_{i}$ of a $d$-dimensional Brownian motion is a one-dimensional Brownian motion, we have the quadratic variation formula $\left[W_{i}, W_{i}\right](t)=t$, which we write informally as

$$
d W_{i}(t) d W_{i}(t)=d t
$$

However, if $i \neq j$, we shall see that independence of $W_{i}$ and $W_{j}$ implies $\left[W_{i}, W_{j}\right](t)=0$, which we write informally as

$$
d W_{i}(t) d W_{j}(t)=0, \quad i \neq j .
$$

We justify this claim.

Let $\Pi=\left\{t_{0}, \ldots, t_{n}\right\}$ be a partition of $[0, T]$. For $i \neq j$, define the sampled cross variation of $W_{i}$ and $W_{j}$ on $[0, T]$ to be

$$
C_{\Pi}=\sum_{k=0}^{n-1}\left[W_{i}\left(t_{k+1}\right)-W_{i}\left(t_{k}\right)\right]\left[W_{j}\left(t_{k+1}\right)-W_{j}\left(t_{k}\right)\right]
$$

The increments appearing on the right-hand side of the equation above are all independent of one another and all have mean zero. Therefore, $\mathbb{E} C_{\Pi}=0$.

We compute $\operatorname{Var}\left(C_{\Pi}\right)$. Note first that

$$
\begin{aligned}
& C_{\Pi}^{2}= \sum_{k=0}^{n-1}\left[W_{i}\left(t_{k+1}\right)-W_{i}\left(t_{k}\right)\right]^{2}\left[W_{j}\left(t_{k+1}\right)-W_{j}\left(t_{k}\right)\right]^{2} \\
&+2 \sum_{\ell<k}^{n-1}\left[W_{i}\left(t_{\ell+1}\right)-W_{i}\left(t_{\ell}\right)\right]\left[W_{\jmath}\left(t_{\ell+1}\right)-W_{j}\left(t_{\ell}\right)\right] \\
& \cdot\left[W_{i}\left(t_{k+1}\right)-W_{i}\left(t_{k}\right)\right]\left[W_{j}\left(t_{k+1}\right)-W_{j}\left(t_{k}\right)\right]
\end{aligned}
$$

All the increments appearing in the sum of cross-terms are independent of one another and all have mean zero. Therefore,

$$
\operatorname{Var}\left(C_{\Pi}\right)=\mathbb{E} C_{\Pi}^{2}=\mathbb{E} \sum_{k=0}^{n-1}\left[W_{i}\left(t_{k+1}\right)-W_{i}\left(t_{k}\right)\right]^{2}\left[W_{j}\left(t_{k+1}\right)-W_{j}\left(t_{k}\right)\right]^{2}
$$

But $\left[W_{i}\left(t_{k+1}\right)-W_{i}\left(t_{k}\right)\right]^{2}$ and $\left[W_{j}\left(t_{k+1}\right)-W_{j}\left(t_{k}\right)\right]^{2}$ are independent of one another, and each has expectation $\left(t_{k+1}-t_{k}\right)$. It follows that

$$
\operatorname{Var}\left(C_{\Pi}\right)=\sum_{k=0}^{n-1}\left(t_{k+1}-t_{k}\right)^{2} \leq\|\Pi\| \cdot \sum_{k=0}^{n-1}\left(t_{k+1}-t_{k}\right)=\|\Pi\| \cdot T .
$$

As $\|\Pi\| \rightarrow 0$, we have $\operatorname{Var}\left(C_{\Pi}\right) \rightarrow 0$, so $C_{\Pi}$ converges to the constant $\mathbb{E} C_{\Pi}=$ 0.

\subsubsection{Itô-Doeblin Formula for Multiple Processes}

To keep the notation as simple as possible, we write the Itô formula for two processes driven by a two-dimensional Brownian motion. In the obvious way, the formula generalizes to any number of processes driven by a Brownian motion of any number (not necessarily the same number) of dimensions. Let $X(t)$ and $Y(t)$ be Itô processes, which means they are processes of the form

$$
\begin{aligned}
& X(t)=X(0)+\int_{0}^{t} \Theta_{1}(u) d u+\int_{0}^{t} \sigma_{11}(u) d W_{1}(u)+\int_{0}^{t} \sigma_{12}(u) d W_{2}(u), \\
& Y(t)=Y(0)+\int_{0}^{t} \Theta_{2}(u) d u+\int_{0}^{t} \sigma_{21}(u) d W_{1}(u)+\int_{0}^{t} \sigma_{22}(u) d W_{2}(u) .
\end{aligned}
$$

The integrands $\Theta_{i}(u)$ and $\sigma_{i j}(u)$ are assumed to be adapted processes. In differential notation, we write

$$
\begin{aligned}
& d X(t)=\Theta_{1}(t) d t+\sigma_{11}(t) d W_{1}(t)+\sigma_{12}(t) d W_{2}(t) \\
& d Y(t)=\Theta_{2}(t) d t+\sigma_{21}(t) d W_{1}(t)+\sigma_{22}(t) d W_{2}(t)
\end{aligned}
$$

The Itô integral $\int_{0}^{t} \sigma_{11}(u) d W_{1}(u)$ accumulates quadratic variation at rate $\sigma_{11}^{2}(t)$ per unit time, and the Itô integral $\int_{0}^{t} \sigma_{12}(u) d W_{2}(u)$ accumulates quadratic variation at rate $\sigma_{12}^{2}(t)$ per unit time. Because both of these integrals appear in $X(t)$, the process $X(t)$ accumulates quadratic variation at rate $\sigma_{11}^{2}(t)+\sigma_{12}^{2}(t)$ per unit time:

$$
[X, X](t)=\int_{0}^{t}\left(\sigma_{11}^{2}(u)+\sigma_{12}^{2}(u)\right) d u .
$$

We may write this equation in differential form as

$$
d X(t) d X(t)=\left(\sigma_{11}^{2}(t)+\sigma_{12}^{2}(t)\right) d t .
$$

One can informally derive (4.6.3) by squaring (4.6.1) and using the multiplication rules

$$
d t d t=0, d t d W_{i}(t)=0, d W_{i}(t) d W_{i}(t)=d t, d W_{i}(t) d W_{j}(t)=0 \text { for } i \neq j .
$$

In a similar way, we may derive the differential formulas

$$
\begin{aligned}
& d Y(t) d Y(t)=\left(\sigma_{21}^{2}(t)+\sigma_{22}^{2}(t)\right) d t \\
& d X(t) d Y(t)=\left(\sigma_{11}(t) \sigma_{21}(t)+\sigma_{12}(t) \sigma_{22}(t)\right) d t
\end{aligned}
$$

Equation (4.6.5) says that, for every $T \geq 0$,

$$
[X, Y](T)=\int_{0}^{T}\left(\sigma_{11}(t) \sigma_{21}(t)+\sigma_{12}(t) \sigma_{22}(t)\right) d t .
$$

The term $[X, Y](T)$ on the left-hand side is defined as follows. Let $\Pi=$ $\left\{t_{0}, t_{1}, \ldots, t_{n}\right\}$ be a partition of $[0, T]$ (i.e., $0=t_{0}<t_{1}<\cdots<t_{n}=T$ ) and set up the sampled cross variation

$$
\sum_{k=0}^{n-1}\left[X\left(t_{k+1}\right)-X\left(t_{k}\right)\right]\left[Y\left(t_{k+1}\right)-Y\left(t_{k}\right)\right]
$$

Now let the number of partition points $n$ go to infinity as the length of the longest subinterval $\|\Pi\|=\max _{0 \leq k \leq n-1}\left(t_{k+1}-t_{k}\right)$ goes to zero. The limit of the sum in (4.6.7) is $[X, Y](T)$. This limit is given by the right-hand side of (4.6.6). The proof of this assertion is similar to the proof of Lemma 4.4.4, with the additional feature that we must use the fact that $\left[W_{1}, W_{2}\right](t)=0$. We omit the details.

The following theorem generalizes the Itô-Doeblin formula of Theorem 4.4.6. The justification, which we omit, is similar to that of Theorem 4.4.6.

Theorem 4.6.2 (Two-dimensional Itô-Doeblin formula). Let $f(t, x, y)$ be a function whose partial derivatives $f_{t}, f_{x}, f_{y}, f_{x x}, f_{x y}, f_{y x}$, and $f_{y y}$ are defined and are continuous. Let $X(t)$ and $Y(t)$ be Itô processes as discussed above. The two-dimensional Itô-Doeblin formula in differential form is

$$
\begin{aligned}
& d f(t, X(t), Y(t)) \\
& =f_{t}(t, X(t), Y(t)) d t+f_{x}(t, X(t), Y(t)) d X(t)+f_{y}(t, X(t), Y(t)) d Y(t) \\
& +\frac{1}{2} f_{x x}(t, X(t), Y(t)) d X(t) d X(t)+f_{x y}(t, X(t), Y(t)) d X(t) d Y(t) \\
& +\frac{1}{2} f_{y y}(t, X(t), Y(t)) d Y(t) d Y(t) \text {. }
\end{aligned}
$$

Before discussing formula (4.6.8), we rewrite it, leaving out $t$ wherever possible, to obtain the same formula in the more compact notation

$$
\begin{array}{rl}
d f(t, X, Y)=f_{t} & d t+f_{x} d X+f_{y} d Y \\
& +\frac{1}{2} f_{x x} d X d X+f_{x y} d X d Y+\frac{1}{2} f_{y y} d Y d Y
\end{array}
$$

The right-hand side of (4.6.9) is the Taylor series expansion of $f$ out to second order. The full expansion would have the additional second-order terms $f_{t t} d t d t, \frac{1}{2} f_{t x} d t d X$, and $\frac{1}{2} f_{t y} d t d Y$, but $d t d t, d t d X$, and $d t d Y$ are zero. The Taylor series expansion actually has two mixed partial terms, $\frac{1}{2} f_{x y} d X d Y$ and $\frac{1}{2} f_{y x} d Y d X$. For functions $f$ whose second partial derivatives exist and are continuous, $f_{x y}=f_{y x}$, and so we have combined these terms into the single term $f_{x y} d X d Y$ in (4.6.9).

The differentials $d X, d Y, d X d X, d X d Y$, and $d Y d Y$ appearing in (4.6.9) are given by (4.6.1)-(4.6.5). Making these substitutions and then integrating (4.6.9), we obtain the Itô-Doeblin formula in integral form:

$$
\begin{aligned}
& f(t, X(t), Y(t))-f(0, X(0), Y(0)) \\
&= \int_{0}^{t}\left[\sigma_{11}(u) f_{x}(u, X(u), Y(u))+\sigma_{21}(u) f_{y}(u, X(u), Y(u))\right] d W_{1}(u) \\
&+\int_{0}^{t}\left[\sigma_{12}(u) f_{x}(u, X(u), Y(u))+\sigma_{22}(u) f_{y}(u, X(u), Y(u))\right] d W_{2}(u)
\end{aligned}
$$



$$
\begin{aligned}
+\int_{0}^{t} & {\left[f_{t}(u, X(u), Y(u))\right.} \\
& +\Theta_{1}(u) f_{x}(u, X(u), Y(u))+\Theta_{2}(u) f_{y}(u, X(u), Y(u)) \\
& +\frac{1}{2}\left(\sigma_{11}^{2}(u)+\sigma_{12}^{2}(u)\right) f_{x x}(u, X(u), Y(u)) \\
& +\left(\sigma_{11}(u) \sigma_{21}(u)+\sigma_{12}(u) \sigma_{22}(u)\right) f_{x y}(u, X(u), Y(u)) \\
& \left.+\frac{1}{2}\left(\sigma_{21}^{2}(u)+\sigma_{22}^{2}(u)\right) f_{y y}(u, X(u), Y(u))\right] d u
\end{aligned}
$$

The right-hand side of this equation has one ordinary (Lebesgue) integral with respect to $d u$ and two Itô integrals, one with respect to $d W_{1}(u)$ and the other with respect to $d W_{2}(u)$. All terms have precise mathematical meanings. This equation demonstrates why it is preferable to work with differential notation, such as in (4.6.9).

Corollary 4.6.3 (Itô product rule). Let $X(t)$ and $Y(t)$ be Itô processes. Then

$$
d(X(t) Y(t))=X(t) d Y(t)+Y(t) d X(t)+d X(t) d Y(t) .
$$

Proof: In (4.6.9), take $f(t, x, y)=x y$, so that $f_{t}=0, f_{x}=y, f_{y}=x$, $f_{x x}=0, f_{x y}=1$, and $f_{y y}=0$.

\subsubsection{Recognizing a Brownian Motion}

A Brownian motion $W(t)$ is a martingale with continuous paths whose quadratic variation is $[W, W](t)=t$. It turns out that these conditions characterize Brownian motion in the sense of the following theorem.

Theorem 4.6.4 (Lévy, one dimension). Let $M(t), t \geq 0$, be a martingale relative to a filtration $\mathcal{F}(t), t \geq 0$. Assume that $M(0)=0, M(t)$ has continuous paths, and $[M, M](t)=t$ for all $t \geq 0$. Then $M(t)$ is a Brownian motion.

IDEA OF THE PROOF: A Brownian motion is a martingale whose increments are normally distributed. The surprising feature of Lévy's Theorem is that the assumptions do not say anything about normality, and yet implicit in the conclusion is the assertion that $M(t)$ is normally distributed.

The method used to establish normality is to first check that in the derivation of the Itô-Doeblin formula, Theorem 4.4.1, for Brownian motion, the only properties of Brownian motion that were used are assumed in this theorem: a continuous process with quadratic variation $[M, M](t)=t$. Therefore, the ItôDoeblin formula may be applied to $M$ with the result that, for any function $f(t, x)$ whose derivatives exist and are continuous,

$$
d f(t, M(t))=f_{t}(t, M(t)) d t+f_{x}(t, M(t)) d M(t)+\frac{1}{2} f_{x x}(t, M(t)) d t .
$$

The last term uses the fact that $d M(t) d M(t)=d t$. In integrated form, (4.6.11) is

$$
\begin{aligned}
f(t, M(t))=f(0, & M(0))+\int_{0}^{t}\left[f_{t}(s, M(s))+\frac{1}{2} f_{x x}(s, M(s))\right] d s \\
& +\int_{0}^{t} f_{x}(s, M(s)) d M(s) .
\end{aligned}
$$

Because $M(t)$ is a martingale, the stochastic integral $\int_{0}^{t} f_{x}(s, M(s)) d M(s)$ is also. (See Exercise $4.1$ for the case of a simple integrand; the general case follows from this exercise upon passage to the limit.) At $t=0$, this stochastic integral takes the value zero, and so its expectation is always zero. Taking expectations in (4.6.12), we obtain

$$
\mathbb{E} f(t, M(t))=f(0, M(0))+\mathbb{E} \int_{0}^{t}\left[f_{t}(s, M(s))+\frac{1}{2} f_{x x}(s, M(s))\right] d s .
$$

We fix a number $u$ and define

$$
f(t, x)=\exp \left\{u x-\frac{1}{2} u^{2} t\right\} .
$$

Then $f_{t}(t, x)=-\frac{1}{2} u^{2} f(t, x), f_{x}(t, x)=u f(t, x)$, and $f_{x x}(t, x)=u^{2} f(t, x)$. In particular,

$$
f_{t}(t, x)+\frac{1}{2} f_{x x}(t, x)=0 .
$$

For this function $f(t, x)$, the second term on the right-hand side of (4.6.13) is zero, and that equation becomes

$$
\mathbb{E} \exp \left\{u M(t)-\frac{1}{2} u^{2} t\right\}=1 .
$$

In other words, we have the moment-generating function formula

$$
\mathbb{E} e^{u M(t)}=e^{\frac{1}{2} u^{2} t} .
$$

This is the moment-generating function for the normal distribution with mean zero and variance $t$ (see (3.2.13)). Hence, that is the distribution that $M(t)$ must have.

The idea used to justify Theorem 4.6.4 can be combined with the twodimensional Itô-Doeblin formula used to show independence. In particular, we have the following two-dimensional version of Lévy's Theorem.

Theorem 4.6.5 (Lévy, two dimensions). Let $M_{1}(t)$ and $M_{2}(t), t \geq 0$, be martingales relative to a filtration $\mathcal{F}(t), t \geq 0$. Assume that for $i=1,2$, we have $M_{i}(0)=0, M_{i}(t)$ has continuous paths, and $\left[M_{i}, M_{i}\right](t)=t$ for all $t \geq 0$. If, in addition, $\left[M_{1}, M_{2}\right](t)=0$ for all $t \geq 0$, then $M_{1}(t)$ and $M_{2}(t)$ are independent Brownian motions. IDEA OF THE PROOF: The one-dimensional Lévy Theorem, Theorem 4.6.4, implies that $M_{1}$ and $M_{2}$ are Brownian motions. To show independence, we examine the joint moment-generating function.

Let $f(t, x, y)$ be a function whose derivatives are defined and continuous. The two-dimensional Itô-Doeblin formula implies that

$$
\begin{aligned}
d f\left(t, M_{1}, M_{2}\right)= & f_{t} d t+f_{x} d M_{1}+f_{y} d M_{2} \\
& \quad+\frac{1}{2} f_{x x} d M_{1} d M_{1}+f_{x y} d M_{1} d M_{2}+f_{y y} d M_{2} d M_{2} \\
= & f_{t} d t+f_{x} d M_{1}+f_{y} d M_{2}+\frac{1}{2} f_{x x} d t+\frac{1}{2} f_{y y} d t,
\end{aligned}
$$

where we have used the assumptions $\left[M_{1}, M_{1}\right](t)=t,\left[M_{2}, M_{2}\right](t)=t$, and $\left[M_{1}, M_{2}\right](t)=$ in their equivalent form $d M_{1}(t) d M_{1}(t)=d t, d M_{2}(t) d M_{2}(t)=$ $d t$, and $d M_{1}(t) d M_{2}(t)=0$. We integrate both sides to obtain

$$
\begin{aligned}
& f\left(t, M_{1}(t), M_{2}(t)\right) \\
& =f\left(0, M_{1}(0), M_{2}(0)\right)+\int_{0}^{t}\left[f_{t}\left(s, M_{1}(s), M_{2}(s)\right)+\frac{1}{2} f_{x x}\left(s, M_{1}(s), M_{2}(s)\right)\right. \\
& \left.+\frac{1}{2} f_{y y}\left(s, M_{1}(s), M_{2}(s)\right)\right] d s \\
& +\int_{0}^{t} f_{x}\left(s, M_{1}(s), M_{2}(s)\right) d M_{1}(s)+\int_{0}^{t} f_{y}\left(x, M_{1}(s), M_{2}(s)\right) d M_{2}(s) .
\end{aligned}
$$

The last two terms on the right-hand side are martingales, starting at zero at time zero, and hence having expectation zero. Therefore,

$$
\begin{array}{r}
\mathbb{E} f\left(t, M_{1}(t), M_{2}(t)\right) \\
=f\left(0, M_{1}(0), M_{2}(0)\right)+\mathbb{E} \int_{0}^{t}\left[f_{t}\left(s, M_{1}(s), M_{2}(s)\right)+\frac{1}{2} f_{x x}\left(s, M_{1}(s), M_{2}(s)\right)\right. \\
\left.\quad+\frac{1}{2} f_{y y}\left(s, M_{1}(s), M_{2}(s)\right)\right] d s .
\end{array}
$$

We now fix numbers $u_{1}$ and $u_{2}$ and define

$$
f(t, x, y)=\exp \left\{u_{1} x+u_{2} y-\frac{1}{2}\left(u_{1}^{2}+u_{2}^{2}\right) t\right\} .
$$

Then $f_{t}(t, x, y)=-\frac{1}{2}\left(u_{1}^{2}+u_{2}^{2}\right) f(t, x, y), f_{x}(t, x, y)=u_{1} f(t, x, y), f_{y}(t, x, y)=$ $u_{2} f(t, x, y) . f_{x x}(t, x, y)=u_{1}^{2} f(t, x, y)$, and $f_{y y}(t, x, y)=u_{2}^{2} f(t, x, y)$. For this function $f(t, x, y)$, the second term on the right-hand side of (4.6.14) is zero. We conclude that

$$
\mathbb{E} \exp \left\{u_{1} M_{1}(t)+u_{2} M_{2}(t)-\frac{1}{2}\left(u_{1}^{2}+u_{2}^{2}\right) t\right\}=1,
$$

which gives us the moment-generating function formula 

$$
\mathbb{E} e^{u_{1} M_{1}(t)+u_{2} M_{2}(t)}=e^{\frac{1}{2} u_{1}^{2} t} \cdot e^{\frac{1}{2} u_{2}^{2} t} .
$$

Because the joint moment-generating function factors into the product of moment-generating functions, $M_{1}(t)$ and $M_{2}(t)$ must be independent.

Example 4.6.6 (Correlated stock prices). Suppose

$$
\begin{aligned}
& \frac{d S_{1}(t)}{S_{1}(t)}=\alpha_{1} d t+\sigma_{1} d W_{1}(t) \\
& \frac{d S_{2}(t)}{S_{2}(t)}=\alpha_{2} d t+\sigma_{2}\left[\rho d W_{1}(t)+\sqrt{1-\rho^{2}} d W_{2}(t)\right],
\end{aligned}
$$

where $W_{1}(t)$ and $W_{2}(t)$ are independent Brownian motions and $\sigma_{1}>0, \sigma_{2}>0$ and $-1 \leq \rho \leq 1$ are constant. To analyze the second stock price process, we define

$$
W_{3}(t)=\rho W_{1}(t)+\sqrt{1-\rho^{2}} W_{2}(t) .
$$

Then $W_{3}(t)$ is a continuous martingale with $W_{3}(0)=0$, and

$$
\begin{aligned}
d W_{3}(t) d W_{3}(t)= & \rho^{2} d W_{1}(t) d W_{1}(t)+2 \rho \sqrt{1-\rho^{2}} d W_{1}(t) d W_{2}(t) \\
& +\left(1-\rho^{2}\right) d W_{2}(t) d W_{2}(t) \\
= & \rho^{2} d t+\left(1-\rho^{2}\right) d t=d t .
\end{aligned}
$$

In other words, $\left[W_{3}, W_{3}\right](t)=t$. According to the one-dimensional Lévy Theorem, Theorem 4.6.4, $W_{3}(t)$ is a Brownian motion. Because we can write the differential of $S_{2}(t)$ as

$$
\frac{d S_{2}(t)}{S_{2}(t)}=\alpha_{2} d t+\sigma_{2} d W_{3}(t),
$$

we see that $S_{2}(t)$ is a geometric Brownian motion with mean rate of return $\alpha_{2}$ and volatility $\sigma_{2}$

The Brownian motions $W_{1}(t)$ and $W_{3}(t)$ are correlated. According to Itô's product rule (Corollary 4.6.3),

$$
\begin{aligned}
d\left(W_{1}(t) W_{3}(t)\right) & =W_{1}(t) d W_{3}(t)+W_{3}(t) d W_{1}(t)+d W_{1}(t) d W_{3}(t) \\
& =W_{1}(t) d W_{3}(t)+W_{3}(t) d W_{1}(t)+\rho d t
\end{aligned}
$$

Integrating, we obtain

$$
W_{1}(t) W_{3}(t)=\int_{0}^{t} W_{1}(s) d W_{3}(s)+\int_{0}^{t} W_{3}(s) d W_{1}(s)+\rho t .
$$

The Itô integrals on the right-hand side have expectation zero, so the covariance of $W_{1}(t)$ and $W_{3}(t)$ is

$$
\mathbb{E}\left[W_{1}(t) W_{3}(t)\right]=\rho t .
$$

Because both $W_{1}(t)$ and $W_{3}(t)$ have standard deviation $\sqrt{t}$, the number $\rho$ is the correlation between $W_{1}(t)$ and $W_{3}(t)$. The case of nonconstant correlation $\rho$ is presented in Exercise $4.17$. 

\subsection{Brownian Bridge}

We conclude this chapter with a the discussion of the Brownian bridge. This is a stochastic process that is like a Brownian motion except that with probability one it reaches a specified point at a specified positive time. We first discuss Gaussian processes in general, the class to which the Brownian bridge belongs, and we then define the Brownian bridge and present its properties. The primary use for the Brownian bridge in finance is as an aid to Monte Carlo simulation. We make no use of it in this text.

\subsubsection{Gaussian Processes}

Definition 4.7.1. A Gaussian process $X(t), t \geq 0$, is a stochastic process that has the property that, for arbitrary times $0<t_{1}<t_{2}<\cdots<t_{n}$, the random variables $X\left(t_{1}\right), X\left(t_{2}\right), \ldots X\left(t_{n}\right)$ are jointly normally distributed.

The joint normal distribution of a set of vectors is determined by their means and covariances. Therefore, for a Gaussian process, the joint distribution of $X\left(t_{1}\right), X\left(t_{2}\right), \ldots, X\left(t_{n}\right)$ is determined by the means and covariances of these random variables. We denote the mean of $X(t)$ by $m(t)$, and, for $s \geq 0$, $t \geq 0$, we denote the covariance of $X(s)$ and $X(t)$ by $c(s, t)$; i.e.,

$$
m(t)=\mathbb{E} X(t), \quad c(s, t)=\mathbb{E}[(X(s)-m(s))(X(t)-m(t))] .
$$

Example 4.7.2 (Brownian motion). Brownian motion $W(t)$ is a Gaussian process. For $0<t_{1}<t_{2}<\cdots<t_{n}$, the increments

$$
I_{1}=W\left(t_{1}\right), I_{2}=W\left(t_{2}\right)-W\left(t_{1}\right), \ldots, I_{n}=W\left(t_{n}\right)-W\left(t_{n-1}\right)
$$

are independent and normally distributed. Writing

$$
W\left(t_{1}\right)=I_{1}, W\left(t_{2}\right)=\sum_{j=1}^{2} I_{j}, \ldots, W\left(t_{n}\right)=\sum_{j=1}^{n} I_{j},
$$

we see that the random variables $W\left(t_{1}\right), W\left(t_{2}\right), \ldots, W\left(t_{n}\right)$ are jointly normally distributed. These random variables are not independent. It is the increments of Brownian motion that are independent. Of course, the mean function for Brownian motion is

$$
m(t)=\mathbb{E} W(t)=0 .
$$

We may compute the covariance by letting $0 \leq s \leq t$ be given and noting that

$$
\begin{aligned}
c(s, t) & =\mathbb{E}[W(s) W(t)] \\
& =\mathbb{E}[W(s)(W(t)-W(s)+W(s))] \\
& =\mathbb{E}[W(s)(W(t)-W(s))]+\mathbb{E}\left[W^{2}(s)\right] .
\end{aligned}
$$

Because $W(s)$ and $W(t)-W(s)$ are independent and both have mean zero, we see that $\mathbb{E}[W(s)(W(t)-W(s))]=0$. The other term, $\mathbb{E}\left[W^{2}(s)\right]$, is the variance of $W(s)$, which is $s$. We conclude that $c(s, t)=s$ when $0 \leq s \leq t$. Reversing the roles of $s$ and $t$, we conclude that $c(s, t)=t$ when $0 \leq t \leq s$. In general, the covariance function for Brownian motion is then

$$
c(s, t)=s \wedge t,
$$

where $s \wedge t$ denotes the minimum of $s$ and $t$.

Example 4.7.3 (Itô integral of a deterministic integrand). Let $\Delta(t)$ be a nonrandom function of time, and define

$$
I(t)=\int_{0}^{t} \Delta(s) d W(s),
$$

where $W(t)$ is a Brownian motion. Then $I(t)$ is a Gaussian process, as we now show.

In the proof of Theorem 4.4.9, we showed that, for fixed $u \in \mathbb{R}$, the process

$$
M_{u}(t)=\exp \left\{u I(t)-\frac{1}{2} u^{2} \int_{0}^{t} \Delta^{2}(s) d s\right\}
$$

is a martingale. We used this fact to argue that

$$
1=M_{u}(0)=\mathbb{E} M_{u}(t)=e^{-\frac{1}{2} u^{2} \int_{0}^{t} \Delta^{2}(s) d s} \cdot \mathbb{E} e^{u I(t)},
$$

and we thus obtained the moment-generating function formula

$$
\mathbb{E} e^{u I(t)}=e^{\frac{1}{2} u^{2} \int_{0}^{t} \Delta^{2}(s) d s} .
$$

The right-hand side is the moment generating function for a normal random variable with mean zero and variance $\int_{0}^{t} \Delta^{2}(s) d s$. Therefore, this is the distribution of $I(t)$.

Although we have shown that $I(t)$ is normally distributed, verification that the process is Gaussian requires more. We must verify that, for $0<t_{1}<$ $t_{2}<\cdots<t_{n}$, the random variables $I\left(t_{1}\right), I\left(t_{2}\right), \ldots, I\left(t_{n}\right)$ are jointly normally distributed. It turns out that the increments

$$
I\left(t_{1}\right)-I(0)=I\left(t_{1}\right), I\left(t_{2}\right)-I\left(t_{1}\right), \ldots, I\left(t_{n}\right)-I\left(t_{n-1}\right)
$$

are normally distributed and independent, and from this the joint normality of $I\left(t_{1}\right), I\left(t_{2}\right), \ldots, I\left(t_{n}\right)$ follows by the same argument as used in Example 4.7.2 for Brownian motion.

We show that, for $0<t_{1}<t_{2}$, the two random increments $I\left(t_{1}\right)-I(0)=$ $I\left(t_{1}\right)$ and $I\left(t_{2}\right)-I\left(t_{1}\right)$ are normally distributed and independent. The argument we provide can be iterated to prove this result for any number of increments. For fixed $u_{2} \in \mathbb{R}$, the martingale property of $M_{u_{2}}$ implies that 

$$
M_{u_{2}}\left(t_{1}\right)=\mathbb{E}\left[M_{u_{2}}\left(t_{2}\right) \mid \mathcal{F}\left(t_{1}\right)\right] .
$$

Now let $u_{1} \in \mathbb{R}$ be fixed. Because $\frac{M_{u_{1}}\left(t_{1}\right)}{M_{u_{2}}\left(t_{1}\right)}$ is $\mathcal{F}\left(t_{1}\right)$-measurable, we may multiply the equation above by this quotient to obtain

$$
\begin{aligned}
& M_{u_{1}}\left(t_{1}\right)= \mathbb{E}\left[\frac{M_{u_{1}}\left(t_{1}\right) M_{u_{2}}\left(t_{2}\right)}{M_{u_{2}}\left(t_{1}\right)} \mid \mathcal{F}\left(t_{1}\right)\right] \\
&= \mathbb{E}\left[\operatorname { e x p } \left\{u_{1} I\left(t_{1}\right)+u_{2}\left(I\left(t_{2}\right)-I\left(t_{1}\right)\right)-\frac{1}{2} u_{1}^{2} \int_{0}^{t_{1}} \Delta^{2}(s) d s\right.\right. \\
&\left.\left.-\frac{1}{2} u_{2}^{2} \int_{t_{1}}^{t_{2}} \Delta^{2}(s) d s\right\} \mid \mathcal{F}\left(t_{1}\right)\right] .
\end{aligned}
$$

We now take expectations

$$
\begin{aligned}
1 & =M_{u_{1}}(0) \\
= & \mathbb{E} M_{u_{1}}\left(t_{1}\right) \\
= & \mathbb{E}\left[\operatorname { e x p } \left\{u_{1} I\left(t_{1}\right)+u_{2}\left(I\left(t_{2}\right)-I\left(t_{1}\right)\right)-\frac{1}{2} u_{1}^{2} \int_{0}^{t_{1}} \Delta^{2}(s) d s\right.\right. \\
& \left.\left.\quad-\frac{1}{2} u_{2}^{2} \int_{t_{1}}^{t_{2}} \Delta^{2}(s) d s\right\}\right] \\
= & \mathbb{E}\left[\exp \left\{u_{1} I\left(t_{1}\right)+u_{2}\left(I\left(t_{2}\right)-I\left(t_{1}\right)\right)\right\}\right] \\
& \cdot \exp \left\{-\frac{1}{2} u_{1}^{2} \int_{0}^{t_{1}} \Delta^{2}(s) d s-\frac{1}{2} u_{2}^{2} \int_{t_{1}}^{t_{2}} \Delta^{2}(s) d s\right\},
\end{aligned}
$$

where we have used the fact that $\Delta^{2}(s)$ is nonrandom to take the integrals of $\Delta^{2}(s)$ outside the expectation on the right-hand side. This leads to the moment-generating function formula

$$
\begin{aligned}
& \mathbb{E}\left[\exp \left\{u_{1} I\left(t_{1}\right)+u_{2}\left(I\left(t_{2}\right)-I\left(t_{1}\right)\right)\right\}\right] \\
& =\exp \left\{\frac{1}{2} u_{1}^{2} \int_{0}^{t_{1}} \Delta^{2}(s) d s\right\} \cdot \exp \left\{\frac{1}{2} u_{2}^{2} \int_{t_{1}}^{t_{2}} \Delta^{2}(s) d s\right\} .
\end{aligned}
$$

The right-hand side is the product of the moment-generating function for a normal random variable with mean zero and variance $\int_{0}^{t_{1}} \Delta^{2}(s) d s$ and the moment-generating function for a normal random variable with mean zero and variance $\int_{t_{1}}^{t_{2}} \Delta^{2}(s) d s$. It follows that $I\left(t_{1}\right)$ and $I\left(t_{2}\right)-I\left(t_{1}\right)$ must have these distributions, and because their joint moment-generating function factors into this product of moment-generating functions, they must be independent.

The covariance of $I\left(t_{1}\right)$ and $I\left(t_{2}\right)$ can be computed using the same trick as in Example 4.7.2 for the covariance of Brownian motion. We have 

$$
\begin{aligned}
c\left(t_{1}, t_{2}\right) & =\mathbb{E}\left[I\left(t_{1}\right) I\left(t_{2}\right)\right] \\
& =\mathbb{E}\left[I\left(t_{1}\right)\left(I\left(t_{2}\right)-I\left(t_{1}\right)+I\left(t_{1}\right)\right)\right] \\
& =\mathbb{E}\left[I\left(t_{1}\right)\left(I\left(t_{2}\right)-I\left(t_{1}\right)\right)\right]+\mathbb{E} I^{2}\left(t_{1}\right) \\
& =\mathbb{E} I\left(t_{1}\right) \cdot \mathbb{E}\left[I\left(t_{2}\right)-I\left(t_{1}\right)\right]+\int_{0}^{t_{1}} \Delta^{2}(s) d s \\
& =\int_{0}^{t_{1}} \Delta^{2}(s) d s .
\end{aligned}
$$

For the general case where $s \geq 0$ and $t \geq 0$ and we do not know the relationship between $s$ and $t$, we have the covariance formula

$$
c(s, t)=\int_{0}^{s \wedge t} \Delta^{2}(u) d u .
$$

\subsubsection{Brownian Bridge as a Gaussian Process}

Definition 4.7.4. Let $W(t)$ be a Brownian motion. Fix $T>0$. We define the Brownian bridge from 0 to 0 on $[0, T]$ to be the process

$$
X(t)=W(t)-\frac{t}{T} W(T), 0 \leq t \leq T .
$$

Note that $\frac{t}{T} W(T)$ as a function of $t$ is the line from $(0,0)$ to $(T, W(T))$. In (4.7.2), we have subtracted this line away from the Brownian motion $W(t)$, so that the resulting process $X(t)$ satisfies

$$
X(0)=X(T)=0 .
$$

Because $W(T)$ enters the definition of $X(t)$ for $0 \leq t \leq T$, the Brownian bridge $X(t)$ is not adapted to the filtration $\mathcal{F}(t)$ generated by $W(t)$. We shall later obtain a different process that has the same distribution as the process $X(t)$ but is adapted to this filtration.

For $0<t_{1}<t_{2}<\cdots<t_{n}<T$, the random variables

$$
X\left(t_{1}\right)=W\left(t_{1}\right)-\frac{t_{1}}{T} W(T), \ldots, X\left(t_{n}\right)=W\left(t_{n}\right)-\frac{t_{n}}{T} W(T)
$$

are jointly normal because $W\left(t_{1}\right), \ldots, W\left(t_{n}\right), W(T)$ are jointly normal. Hence, the Brownian bridge from 0 to 0 is a Gaussian process. Its mean function is easily seen to be

$$
m(t)=\mathbb{E} X(t)=\mathbb{E}\left[W(t)-\frac{t}{T} W(T)\right]=0 .
$$

For $s, t \in(0, T)$, we compute the covariance function 

$$
\begin{aligned}
& c(s, t) \\
& =\mathbb{E}\left[\left(W(s)-\frac{s}{T} W(T)\right)\left(W(t)-\frac{t}{T} W(T)\right)\right] \\
& =\mathbb{E}[W(s) W(t)]-\frac{t}{T} \mathbb{E}[W(s) W(T)]-\frac{s}{T} \mathbb{E}[W(t) W(T)]+\frac{s t}{T^{2}} \mathbb{E} W^{2}(T) \\
& =s \wedge t-\frac{2 s t}{T}+\frac{s t}{T}=s \wedge t-\frac{s t}{T} .
\end{aligned}
$$

Definition 4.7.5. Let $W(t)$ be a Brownian motion. Fix $T>0, a \in \mathbb{R}$, and $b \in \mathbb{R}$. We define the Brownian bridge from $a$ to $b$ on $[0, T]$ to be the process

$$
X^{a \rightarrow b}(t)=a+\frac{(b-a) t}{T}+X(t), 0 \leq t \leq T,
$$

where $X(t)=X^{0 \rightarrow 0}$ is the Brownian bridge from 0 to 0 of Definition 4.7.4.

The function $a+\frac{(b-a) t}{T}$, as a function of $t$, is the line from $(0, a)$ to $(T, b)$. When we add this line to the Brownian bridge from 0 to 0 on $[0, T]$, we obtain a process that begins at $a$ at time 0 and ends at $b$ at time $T$. Adding a nonrandom function to a Gaussian process gives us another Gaussian process. The mean function is affected:

$$
m^{a \rightarrow b}(t)=\mathbb{E} X^{a \rightarrow b}(t)=a+\frac{(b-a) t}{T} .
$$

However, the covariance function is not affected:

$$
c^{a \rightarrow b}(s, t)=\mathbb{E}\left[\left(X^{a \rightarrow b}(s)-m^{a \rightarrow b}(s)\right)\left(X^{a \rightarrow b}(t)-m^{a \rightarrow b}(t)\right)\right]=s \wedge t-\frac{s t}{T} .
$$

\subsubsection{Brownian Bridge as a Scaled Stochastic Integral}

We cannot write the Brownian bridge as a stochastic integral of a deterministic integrand because the variance of the Brownian bridge,

$$
\mathbb{E} X^{2}(t)=c(t, t)=t-\frac{t^{2}}{T}=\frac{t(T-t)}{T},
$$

increases for $0 \leq t \leq \frac{T}{2}$ and then decreases for $\frac{T}{2} \leq t \leq T$. In Example 4.7.3, the variance of $I(t)=\int_{0}^{t} \Delta(u) d W(u)$ is $\int_{0}^{t} \Delta^{2}(u) d u$, which is nondecreasing in $t$. However, we can obtain a process with the same distribution as the Brownian bridge from 0 to 0 as a scaled stochastic integral. In particular, consider

$$
Y(t)=(T-t) \int_{0}^{t} \frac{1}{T-u} d W(u), 0 \leq t<T .
$$

The integral

$$
I(t)=\int_{0}^{t} \frac{1}{T-u} d W(u)
$$

is a Gaussian process of the type discussed in Example 4.7.3, provided $t<T$ so the integrand is defined. For $0<t_{1}<t_{2}<\cdots<t_{n}<T$, the random variables

$$
Y\left(t_{1}\right)=\left(T-t_{1}\right) I\left(t_{1}\right), Y\left(t_{2}\right)=\left(T-t_{2}\right) I\left(t_{2}\right), \ldots, Y\left(t_{n}\right)=\left(T-t_{n}\right) I\left(t_{n}\right)
$$

are jointly normal because $I\left(t_{1}\right), I\left(t_{2}\right), \ldots, I\left(t_{n}\right)$ are jointly normal. In particular, $Y$ is a Gaussian process.

The mean and covariance functions of $I$ are

$$
\begin{aligned}
m^{I}(t) & =0 \\
c^{I}(s, t) & =\int_{0}^{s \wedge t} \frac{1}{(T-u)^{2}} d u=\frac{1}{T-s \wedge t}-\frac{1}{T} \text { for all } s, t \in[0, T) .
\end{aligned}
$$

This means that the mean function for $Y$ is $m^{Y}(t)=0$. To compute the covariance function for $Y$, we assume for the moment that $0 \leq s \leq t<T$ so that

$$
c^{I}(s, t)=\frac{1}{T-s}-\frac{1}{T}=\frac{s}{T(T-s)} .
$$

Then

$$
\begin{aligned}
c^{Y}(s, t) & =\mathbb{E}[(T-s)(T-t) I(s) I(t)] \\
& =(T-s)(T-t) \frac{s}{T(T-s)} \\
& =\frac{(T-t) s}{T} \\
& =s-\frac{s t}{T} .
\end{aligned}
$$

If we had taken $0 \leq t \leq s<T$, the roles of $s$ and $t$ would have been reversed. In general,

$$
c^{Y}(s, t)=s \wedge t-\frac{s t}{T} \text { for all } s, t \in[0, T) .
$$

This is the same covariance formula (4.7.3) we obtained for the Brownian bridge. Because the mean and covariance functions for a Gaussian process completely determine the distribution of the process, we conclude that the process $Y$ has the same distribution as the Brownian bridge from 0 to 0 on $[0, T]$

We now consider the variance

$$
\mathbb{E} Y^{2}(t)=c^{Y}(t, t)=\frac{t(T-t)}{T}, 0<t<T .
$$

Note that, as $t \uparrow T$, this variance converges to 0 . In other words, as $t \uparrow$ $T$, the random process $Y(t)$, which always has mean zero, has a variance that converges to zero. We did not initially define $Y(T)$, but this observation suggests that it makes sense to define $Y(T)=0$. If we do that, then $Y(t)$ is continuous at $t=T$. We summarize this discussion with the following theorem. Theorem 4.7.6. Define the process

$$
Y(t)= \begin{cases}(T-t) \int_{0}^{t} \frac{1}{T-u} d W(u) & \text { for } 0 \leq t<T, \\ 0 & \text { for } t=T .\end{cases}
$$

Then $Y(t)$ is a continuous Gaussian process on $[0, T]$ and has mean and covariance functions

$$
\begin{aligned}
m^{Y}(t) & =0, t \in[0, T] \\
c^{Y}(s, t) & =s \wedge t-\frac{s t}{T} \text { for all } s, t \in[0, T] .
\end{aligned}
$$

In particular, the process $Y(t)$ has the same distribution as the Brownian bridge from 0 to 0 on $[0, T]$ (Definition 4.7.5).

We note that the process $Y(t)$ is adapted to the filtration generated by the Brownian motion $W(t)$. It is interesting to compute the stochastic differential of $Y(t)$, which is

$$
\begin{aligned}
d Y(t) & =\int_{0}^{t} \frac{1}{T-u} d W(u) \cdot d(T-t)+(T-t) \cdot d \int_{0}^{t} \frac{1}{T-u} d W(u) \\
& =-\int_{0}^{t} \frac{1}{T-u} d W(u) \cdot d t+d W(t) \\
& =-\frac{Y(t)}{T-t} d t+d W(t) .
\end{aligned}
$$

If $Y(t)$ is positive as $t$ approaches $T$, the drift term $-\frac{Y(t)}{T-t} d t$ becomes large in absolute value and is negative. This drives $Y(t)$ toward zero. On the other hand, if $Y(t)$ is negative, the drift term becomes large and positive, and this again drives $Y(t)$ toward zero. This strongly suggests, and it is indeed true, that as $t \uparrow T$ the process $Y(t)$ converges to zero almost surely.

\subsubsection{Multidimensional Distribution of the Brownian Bridge}

We fix $a \in \mathbb{R}$ and $b \in \mathbb{R}$ and let $X^{a \rightarrow b}(t)$ denote the Brownian bridge from $a$ to $b$ on $[0, T]$. We also fix $0=t_{0}<t_{1}<t_{2}<\cdots<t_{n}<T$. In this section, we compute the joint density of $X^{a \rightarrow b}\left(t_{1}\right), \ldots, X^{a \rightarrow b}\left(t_{n}\right)$.

We recall that the Brownian bridge from $a$ to $b$ has the mean function

$$
m^{a \rightarrow b}(t)=a+\frac{(b-a) t}{T}=\frac{(T-t) a}{T}+\frac{b t}{T}
$$

and covariance function

$$
c(s, t)=s \wedge t-\frac{s t}{T} .
$$

When $s \leq t$, we may write this as 

$$
c(s, t)=s-\frac{s t}{T}=\frac{s(T-t)}{T}, 0 \leq s \leq t \leq T .
$$

To simplify notation, we set $\tau_{j}=T-t_{j}$ so that $\tau_{0}=T$. We define random variables

$$
Z_{j}=\frac{X^{a \rightarrow b}\left(t_{j}\right)}{\tau_{j}}-\frac{X^{a \rightarrow b}\left(t_{j-1}\right)}{\tau_{j-1}} .
$$

Because $X^{a \rightarrow b}\left(t_{1}\right), \ldots, X^{a \rightarrow b}\left(t_{n}\right)$ are jointly normal, so are $Z\left(t_{1}\right), \ldots, Z\left(t_{n}\right)$. We compute

$$
\begin{aligned}
\mathbb{E} Z_{j} & =\frac{1}{\tau_{j}} \mathbb{E} X^{a \rightarrow b}\left(t_{j}\right)-\frac{1}{\tau_{j}} \mathbb{E} X^{a \rightarrow b}\left(t_{j+1}\right) \\
& =\frac{a}{T}+\frac{b t_{j}}{T \tau_{j}}-\frac{a}{T}-\frac{b t_{j-1}}{T \tau_{j-1}} \\
& =\frac{b t_{j}\left(T-t_{j-1}\right)-b t_{j-1}\left(T-t_{j}\right)}{T \tau_{j} \tau_{j-1}} \\
& =\frac{b\left(t_{j}-t_{j-1}\right)}{\tau_{j} \tau_{j-1}}
\end{aligned}
$$

Furthermore,

$$
\begin{aligned}
\operatorname{Var}\left(Z_{j}\right) & =\frac{1}{\tau_{j}^{2}} \operatorname{Var}\left(X^{a \rightarrow b}\left(t_{j}\right)\right)-\frac{2}{\tau_{j} \tau_{j-1}} \operatorname{Cov}\left(X^{a \rightarrow b}\left(t_{j}\right), X^{a \rightarrow b}\left(t_{j-1}\right)\right) \\
& \quad+\frac{1}{\tau_{j-1}^{2}} \operatorname{Var}\left(X^{a \rightarrow b}\left(t_{j-1}\right)\right) \\
& =\frac{1}{\tau_{j}^{2}} c\left(t_{j}, t_{j}\right)-\frac{2}{\tau_{j} \tau_{j-1}} c\left(t_{j}, t_{j-1}\right)+\frac{1}{\tau_{j-1}^{2}} c\left(t_{j-1}, t_{j-1}\right) \\
& =\frac{t_{j}}{T \tau_{j}}-\frac{2 t_{j-1}}{T \tau_{j-1}}+\frac{t_{j-1}}{T \tau_{j-1}} \\
& =\frac{t_{j}\left(T-t_{j-1}\right)-2 t_{j-1}\left(T-t_{j}\right)+t_{j-1}\left(T-t_{j}\right)}{T \tau_{j} \tau_{j-1}} \\
& =\frac{t_{j}-t_{j-1}}{\tau_{j} \tau_{j-1}} .
\end{aligned}
$$

Finally, we compute the covariance of $Z_{i}$ and $Z_{j}$ when $i<j$. We obtain

$$
\begin{aligned}
\operatorname{Cov}\left(Z_{i}, Z_{j}\right)= & \frac{1}{\tau_{i} \tau_{j}} c\left(t_{i}, t_{j}\right)-\frac{1}{\tau_{i} \tau_{j-1}} c\left(t_{i}, t_{j-1}\right)-\frac{1}{\tau_{i-1} \tau_{j}} c\left(t_{i-1}, t_{j}\right) \\
& +\frac{1}{\tau_{i-1} \tau_{j-1}} c\left(t_{i-1}, t_{j-1}\right) \\
= & \frac{t_{i}\left(T-t_{j}\right)}{T \tau_{i} \tau_{j}}-\frac{t_{i}\left(T-t_{j-1}\right)}{T \tau_{i} \tau_{j-1}}-\frac{t_{i-1}\left(T-t_{j}\right)}{T \tau_{i-1} \tau_{j}}+\frac{t_{i-1}\left(T-t_{j-1}\right)}{T \tau_{i-1} \tau_{j-1}} \\
= & 0
\end{aligned}
$$

We conclude that the normal random variables $Z_{1}, \ldots, Z_{n}$ are independent, and we can write down their joint density, which is

$$
\begin{aligned}
& f_{Z\left(t_{1}\right), \ldots, Z\left(t_{n}\right)}\left(z_{1}, \ldots, z_{n}\right) \\
&=\prod_{j=1}^{n} \frac{1}{\sqrt{2 \pi \frac{t_{j}-t_{j-1}}{\tau_{j} \tau_{j-1}}}} \exp \left\{-\frac{1}{2} \cdot \frac{\left(z_{j}-\frac{b\left(t_{j}-t_{j-1}\right)}{\tau_{j} \tau_{j-1}}\right)^{2}}{\frac{t_{j}-t_{j-1}}{\tau_{j} \tau_{j-1}}}\right\} \\
&=\exp \left\{-\frac{1}{2} \sum_{j=1}^{n} \frac{\left(z_{j}-\frac{b\left(t_{j}-t_{j-1}\right)}{\tau_{j} \tau_{j-1}}\right)^{2}}{\frac{t_{j}-t_{j-1}}{\tau_{j} \tau_{j-1}}}\right\} \cdot \prod_{j=1}^{n} \frac{1}{\sqrt{2 \pi \frac{t_{j}-t_{j-1}}{\tau_{j} \tau_{j-1}}}}
\end{aligned}
$$

We make the change of variables

$$
z_{j}=\frac{x_{j}}{\tau_{j}}-\frac{x_{j-1}}{\tau_{j-1}}, \quad j=1, \ldots, n,
$$

where $x_{0}=a$, to find the joint density for $X^{a \rightarrow b}\left(t_{1}\right), \ldots, X^{a \rightarrow b}\left(t_{n}\right)$. We work first on the sum in the exponent to see the effect of this change of variables. We have

$$
\begin{aligned}
& \sum_{j=1}^{n} \frac{\left(z_{j}-\frac{b\left(t_{j}-t_{j-1}\right)}{\tau_{j} \tau_{j-1}}\right)^{2}}{\frac{t_{j}-t_{j-1}}{\tau_{j} \tau_{j-1}}} \\
& =\sum_{j=1}^{n} \frac{\tau_{j} \tau_{j-1}}{t_{j}-t_{j-1}}\left(\frac{x_{j}}{\tau_{j}}-\frac{x_{j-1}}{\tau_{j-1}}-\frac{b\left(t_{j}-t_{j-1}\right.}{\tau_{j} \tau_{j-1}}\right)^{2} \\
& =\sum_{j=1}^{n} \frac{\tau_{j} \tau_{j-1}}{t_{j}-t_{j-1}}\left(\frac{x_{j}^{2}}{\tau_{j}^{2}}+\frac{x_{j-1}^{2}}{\tau_{j-1}^{2}}+\frac{b^{2}\left(t_{j}-t_{j-1}\right)^{2}}{\tau_{j}^{2} \tau_{j-1}^{2}}-\frac{2 x_{j} x_{j-1}}{\tau_{j} \tau_{j-1}}\right. \\
& \left.\quad-\frac{2 x_{j} b\left(t_{j}-t_{j-1}\right)}{\tau_{j}^{2} \tau_{j-1}}+\frac{2 x_{j-1} b\left(t_{j}-t_{j-1}\right)}{\tau_{j} \tau_{j-1}^{2}}\right) \\
& =\sum_{j=1}^{n}\left(\frac{\tau_{j-1} x_{j}^{2}}{\tau_{j}\left(t_{j}-t_{j-1}\right)}+\frac{\tau_{j} x_{j-1}^{2}}{\tau_{j-1}\left(t_{j}-t_{j-1}\right)}+\frac{b^{2}\left(t_{j}-t_{j-1}\right)}{\tau_{j} \tau_{j-1}}-\frac{2 x_{j} x_{j-1}}{t_{j}-t_{j-1}}\right. \\
& \left.\tau_{j}+\frac{2 x_{j-1} b}{\tau_{j-1}}\right) \\
& =\sum_{j=1}^{n}\left[\frac{x_{j}^{2}}{t_{j}-t_{j-1}}\left(1+\frac{\tau_{j-1}-\tau_{j}}{\tau_{j}}\right)+\frac{x_{j-1}^{2}}{t_{j}-t_{j-1}}\left(1-\frac{\tau_{j-1}-\tau_{j}}{\tau_{j-1}}\right)\right. \\
& \left.-\frac{2 x_{j} x_{j-1}}{t_{j}-t_{j-1}}\right]+b^{2} \sum_{j=1}^{n}\left(\frac{1}{\tau_{j}}-\frac{1}{\tau_{j-1}}\right)-2 b \sum_{j=1}^{n}\left(\frac{x_{j}}{\tau_{j}}-\frac{x_{j-1}}{\tau_{j-1}}\right)
\end{aligned}
$$

Now

$$
\tau_{j-1}-\tau_{j}=\left(T-t_{j-1}\right)-\left(T-t_{j}\right)=t_{j}-t_{j-1},
$$

and so this last expression is equal to

$$
\begin{gathered}
\sum_{j=1}^{n}\left[\frac{\left.x_{j}^{2}-2 x_{j} x_{j-1}+x_{j-1}^{2}\right]}{t_{j}-t_{j-1}}+\sum_{j=1}^{n}\left(\frac{x_{j}^{2}}{\tau_{j}}-\frac{x_{j-1}^{2}}{\tau_{j-1}}\right)\right. \\
+b^{2} \sum_{j=1}^{n}\left(\frac{1}{\tau_{j}}-\frac{1}{\tau_{j-1}}\right)-2 b \sum_{j=1}^{n}\left(\frac{x_{j}}{\tau_{j}}-\frac{x_{j-1}}{\tau_{j-1}}\right) \\
=\sum_{j=1}^{n} \frac{\left(x_{j}-x_{j-1}\right)^{2}}{t_{j}-t_{j-1}}+\frac{x_{n}^{2}}{T-t_{n}}-\frac{a^{2}}{T}+b^{2}\left(\frac{1}{T-t_{n}}-\frac{1}{T}\right) \\
-2 b\left(\frac{x_{n}}{T-t_{n}}-\frac{a}{T}\right) \\
=\sum_{j=1}^{n} \frac{\left(x_{j}-x_{j-1}\right)^{2}}{t_{j}-t_{j-1}}+\frac{\left(b-x_{n}\right)^{2}}{T-t_{n}}-\frac{(b-a)^{2}}{T} .
\end{gathered}
$$

In conclusion, when we change variables from $z_{j}$ to $x_{j}$, we have the equation

$$
\begin{aligned}
& \exp \left\{-\frac{1}{2} \sum_{j=1}^{n} \frac{\left(z_{j}-\frac{b\left(t_{j}-t_{j}-1\right)}{\tau_{j} \tau_{j-1}}\right)^{2}}{\frac{t_{j}-t_{j-1}}{\tau_{j} \tau_{j-1}}}\right\} \\
&=\exp \left\{-\frac{1}{2} \sum_{j=1}^{n} \frac{\left(x_{j}-x_{j-1}\right)^{2}}{t_{j}-t_{j-1}}-\frac{\left(b-x_{n}\right)^{2}}{2\left(T-t_{n}\right)}+\frac{(b-a)^{2}}{2 T}\right\} .
\end{aligned}
$$

To change a density, we also need to account for the Jacobian of the change of variables. In this case, we have

$$
\begin{aligned}
\frac{\partial z_{j}}{\partial x_{j}} & =\frac{1}{\tau_{j}}, j=1, \ldots, n, \\
\frac{\partial z_{j}}{\partial x_{j-1}} & =-\frac{1}{\tau_{j-1}}, j=2, \ldots, n,
\end{aligned}
$$

and all other partial derivatives are zero. This leads to the Jacobian matrix

$$
J=\left[\begin{array}{cccc}
\frac{1}{\tau_{1}} & 0 & \cdots & 0 \\
-\frac{1}{\tau_{1}} & \frac{1}{\tau_{2}} & \cdots & 0 \\
\vdots & \vdots & & \vdots \\
0 & 0 & \cdots & \frac{1}{\tau_{n}}
\end{array}\right],
$$

whose determinant is $\prod_{j=1}^{n} \frac{1}{\tau_{j}}$. Multiplying $f_{Z\left(t_{1}\right), \ldots, Z\left(t_{n}\right)}\left(z_{1}, \ldots, z_{n}\right)$ by this determinant and using the change of variables worked out above, we obtain the density for $X^{a \rightarrow b}\left(t_{1}\right), \ldots, X^{a \rightarrow b}\left(t_{n}\right)$, 

$$
\begin{aligned}
& f_{X^{a \rightarrow b}\left(t_{1}\right), \ldots, X^{a \rightarrow b}\left(t_{n}\right)}\left(x_{1}, \ldots, x_{n}\right) \\
& =\prod_{j=1}^{n} \frac{1}{\sqrt{2 \pi\left(t_{j}-t_{j-1}\right)}} \sqrt{\frac{\tau_{j-1}}{\tau_{j}}} \\
& \cdot \exp \left\{-\frac{1}{2} \sum_{j=1}^{n} \frac{\left(x_{j}-x_{j-1}\right)^{2}}{t_{j}-t_{j-1}}-\frac{\left(b-x_{n}\right)^{2}}{2\left(T-t_{n}\right)}+\frac{(b-a)^{2}}{2 T}\right\} \\
& =\sqrt{\frac{T}{T-t_{n}}} \cdot \prod_{j=1}^{n} \frac{1}{\sqrt{2 \pi\left(t_{j}-t_{j-1}\right)}} \\
& \cdot \exp \left\{-\frac{1}{2} \sum_{j=1}^{n} \frac{\left(x_{\jmath}-x_{j-1}\right)^{2}}{t_{j}-t_{j-1}}-\frac{\left(b-x_{n}\right)^{2}}{2\left(T-t_{n}\right)}+\frac{(b-a)^{2}}{2 T}\right\} \\
& =\frac{p\left(T-t_{n}, x_{n}, b\right)}{p(T, a, b)} \prod_{j=1}^{n} p\left(t_{j}-t_{j-1}, x_{j-1}, x_{j}\right) \text {, }
\end{aligned}
$$

where

$$
p(\tau, x, y)=\frac{1}{\sqrt{2 \pi \tau}} \exp \left\{-\frac{(y-x)^{2}}{2 \tau}\right\}
$$

is the transition density for Brownian motion.

\subsubsection{Brownian Bridge as a Conditioned Brownian Motion}

The joint density (4.7.6) for $X^{a \rightarrow b}\left(t_{1}\right), \ldots, X^{a \rightarrow b}\left(t_{n}\right)$ permits us to give one more interpretation for the Brownian bridge from $a$ to $b$ on $[0, T]$. It is a Brownian motion $W(t)$ on this time interval, starting at $W(0)=a$ and conditioned to arrive at $b$ at time $T$ (i.e., conditioned on $W(T)=b$ ). Let $0=t_{0}<t_{1}<t_{2}<\cdots<t_{n}<T$ be given. The joint density of $W\left(t_{1}\right), \ldots, W\left(t_{n}\right), W(T)$ is

$$
f_{W\left(t_{1}\right), \ldots, W\left(t_{n}\right), W(T)}\left(x_{1}, \ldots, x_{n}, b\right)=p\left(T-t_{n}, x_{n}, b\right) \prod_{j=1}^{n} p\left(t_{j}-t_{j-1}, x_{j-1}, x_{j}\right)
$$

where $W(0)=x_{0}=a$. This is because $p\left(t_{1}-t_{0}, x_{0}, x_{1}\right)=p\left(t_{1}, a, x_{1}\right)$ is the density for the Brownian motion going from $W(0)=a$ to $W\left(t_{1}\right)=x_{1}$ in the time between $t=0$ and $t=t_{1}$. Similarly, $p\left(t_{2}-t_{1}, x_{1}, x_{2}\right)$ is the density for going from $W\left(t_{1}\right)=x_{1}$ to $W\left(t_{2}\right)=x_{2}$ between time $t=t_{1}$ and $t=t_{2}$. The joint density for $W\left(t_{1}\right)$ and $W\left(t_{2}\right)$ is then the product

$$
p\left(t_{1}, a, x_{1}\right) p\left(t_{2}-t_{1}, x_{1}, x_{2}\right) .
$$

Continuing in this way, we obtain the joint density (4.7.7). The marginal density of $W(T)$ is $p(T, a, b)$. The density of $W\left(t_{1}\right), \ldots, W\left(t_{n}\right)$ conditioned on $W(T)=b$ is thus the quotient 

$$
\frac{p\left(T-t_{n}, x_{n}, b\right)}{p(T, a, b)} \prod_{j=1}^{n} p\left(t_{j}-t_{j-1}, x_{j-1}, x_{j}\right)
$$

and this is $f_{X^{a \rightarrow b}\left(t_{1}\right), \ldots, X^{a \rightarrow b}\left(t_{n}\right)}\left(x_{1}, \ldots, x_{n}\right)$ of (4.7.6).

Finally, let us define

$$
M^{a \rightarrow b}(T)=\max _{0 \leq t \leq T} X^{a \rightarrow b}(t)
$$

to be the maximum value obtained by the Brownian bridge from $a$ to $b$ on $[0, T]$. This random variable has the following distribution.

Corollary 4.7.7. The density of $M^{a \rightarrow b}(T)$ is

$$
f_{M^{a \rightarrow b}(T)}(y)=\frac{2(2 y-b-a)}{T} e^{-\frac{2}{T}(y-a)(y-b)}, y>\max \{a, b\} .
$$

Proof: Because the Brownian bridge from 0 to $w$ on $[0, T]$ is a Brownian motion conditioned on $W(T)=w$, the maximum of $X^{0 \rightarrow w}$ on $[0, T]$ is the maximum of $W$ on $[0, T]$ conditioned on $W(T)=w$. Therefore, the density of $M^{0 \rightarrow w}(T)$ was computed in Corollary 3.7.4 and is

$$
f_{M^{0 \rightarrow w}(T)}(m)=\frac{2(2 m-w)}{T} e^{-\frac{2 m(m-w)}{T}}, w<m, m>0 .
$$

The density of $f_{M^{a \rightarrow b}(T)}(y)$ can be obtained by translating from the initial condition $W(0)=a$ to $W(0)=0$ and using (4.7.9). In particular, in (4.7.9) we replace $m$ by $y-a$ and replace $w$ by $b-a$. This results in (4.7.8).

\subsection{Summary}

Let $W(t)$ be a Brownian motion and $\Delta(t)$ a stochastic process adapted to the filtration of the Brownian motion. The Itô integral

$$
I(t)=\int_{0}^{t} \Delta(u) d W(u)
$$

is a martingale. Because it is zero at time $t=0$, its expectation is zero for all $t$. Its variance is given by Itô's isometry

$$
\mathbb{E} I^{2}(t)=\mathbb{E} \int_{0}^{t} \Delta^{2}(u) d u .
$$

The quadratic variation accumulated by the Itô integral up to time $t$ is

$$
[I, I](t)=\int_{0}^{t} \Delta^{2}(u) d u .
$$

These assertions appear in Theorem 4.3.1. Note that the quadratic variation (4.8.3) is computed path-by-path and the result may depend on the path, whereas the variance (4.8.2) is an average over all paths. In differential notation, we write (4.8.1) as

$$
d I(t)=\Delta(t) d W(t)
$$

and $(4.8 .3)$ as

$$
d I(t) d I(t)=\Delta^{2}(t) d W(t) d W(t)=\Delta^{2}(t) d t .
$$

An Itô process (Definition 4.4.3) is a process of the form

$$
X(t)=X(0)+\int_{0}^{t} \Delta(u) d W(u)+\int_{0}^{t} \Theta(u) d u
$$

where $X(0)$ is nonrandom and $\Delta(u)$ and $\Theta(u)$ are adapted stochastic processes. According to Lemma 4.4.4, the quadratic variation accumulated by $X$ up to time $t$ is

$$
[X, X](t)=\int_{0}^{t} \Delta^{2}(u) d u .
$$

In differential notation, we write (4.8.4) as

$$
d X(t)=\Delta(t) d W(t)+\Theta(t) d t
$$

and (4.8.5) as

$$
\begin{aligned}
d X(t) d X(t) & =(\Delta(t) d W(t)+\Theta(t) d t)^{2} \\
& =\Delta^{2}(t) d W(t) d W(t)+2 \Delta(t) \Theta(t) d W(t) d t+\Theta^{2}(t) d t d t \\
& =\Delta^{2}(t) d t
\end{aligned}
$$

where we have used the multiplication table

$$
d W(t) d W(t)=d t, \quad d W(t) d t=d t d W(t)=0, \quad d t d t=0 .
$$

Suppose $X$ and $Y$ are Itô processes with differentials

$$
\begin{aligned}
& d X(t)=\Theta_{1}(t) d t+\sigma_{11}(t) d W_{1}(t)+\sigma_{12}(t) d W_{2}(t) \\
& d Y(t)=\Theta_{2}(t) d t+\sigma_{21}(t) d W_{1}(t)+\sigma_{22}(t) d W_{2}(t)
\end{aligned}
$$

where $W_{1}$ and $W_{2}$ are independent Brownian motions. Then

$$
\begin{aligned}
& d X(t) d X(t)=\left(\sigma_{11}^{2}(t)+\sigma_{12}^{2}(t)\right) d t \\
& d X(t) d Y(t)=\left(\sigma_{11}(t) \sigma_{21}(t)+\sigma_{12}(t) \sigma_{22}(t)\right) d t \\
& d Y(t) d Y(t)=\left(\sigma_{21}^{2}(t)+\sigma_{22}^{2}(t)\right) d t .
\end{aligned}
$$

Equations (4.8.8)-(4.8.10) can be obtained by multiplying the equations (4.8.6) and (4.8.7) for $d X(t)$ and $d Y(t)$ and using the multiplication table 

$$
d W_{i}(t) d W_{i}(t)=d t, d W_{i}(t) d t=d t d W_{i}(t)=0, d t d t=0
$$

and

$$
d W_{1}(t) d W_{2}(t)=0 .
$$

Equation (4.8.11) holds for independent Brownian motions. If instead we had

$$
d W_{1}(t) d W_{2}(t)=\rho d t
$$

for a constant $\rho \in[-1,1]$, then $\rho$ would be the correlation between $W_{1}(t)$ and $W_{2}(t)$ (i.e., $\left.\mathbb{E}\left[W_{1}(t) W_{2}(t)\right]=\rho t\right)$.

Now suppose $f(t, x, y)$ is a function of the time variable $t$ and two dummy variables $x$ and $y$. The multidimensional Itô-Doeblin formula (Theorem 4.6.2) says

$$
\begin{aligned}
& d f(t, X(t), Y(t)) \\
& \begin{aligned}
=f_{t}(t, & X(t), Y(t)) d t+f_{x}(t, X(t), Y(t)) d X(t)+f_{y}(t, X(t), Y(t)) d Y(t) \\
& \frac{1}{2} f_{x x}(t, X(t), Y(t)) d X(t) d X(t)+f_{x y}(t, X(t), Y(t)) d X(t) d Y(t) \\
& +\frac{1}{2} f_{y y}(t, X(t), Y(t)) d Y(t) d Y(t) .
\end{aligned}
\end{aligned}
$$

Replacing all the differentials on the right-hand side of (4.8.12) by their formulas (4.8.6)-(4.8.10) and integrating, one obtains a formula for the stochastic process $f(t, X(t), Y(t))$ as the sum of $f(0, X(0), Y(0))$, an ordinary integral with respect to time, an Itô integral with respect to $d W_{1}$, and an Itô integral with respect to $d W_{2}$.

There are two important special cases of (4.8.12). If the second process $Y$ is not present, (4.8.12) reduces to the Itô-Doeblin formula for one process (Theorem 4.4.6):

$d f(t, X(t))=f_{t}(t, X(t)) d t+f_{x}(t, X(t)) d X(t)+\frac{1}{2} f_{x x}(t, X(t)) d X(t) d X(t)$.

If both $X$ and $Y$ are present and $f(t, x, y)=x y$, then (4.8.12) gives us Itô's product rule (Corollary 4.6.3):

$$
d(X(t) Y(t))=X(t) d Y(t)+Y(t) d X(t)+d X(t) d Y(t) .
$$

Using the Itô-Doeblin formula, we can derive the Black-Scholes-Merton partial differential equation. This was done in Section 4.5, and that section is summarized here. Let the stock price $S(t)$ be a geometric Brownian motion:

$$
d S(t)=\alpha S(t) d t+\sigma S(t) d W(t) .
$$

Let $c(t, S(t))$ be the price at time $t \in[0, T]$ of a European call paying $(S(T)-$ $K)^{+}$at expiration time $T$. Suppose we sell this call for $X(0)=c(0, S(0))$ at time zero and, starting with initial capital $X(0)$, invest in a stock and a money market account paying a constant rate of interest $r$. If $\Delta(t)$ is the number of shares of stock held by the portfolio at time $t$, then

$$
d X(t)=\Delta(t) d S(t)+r(X(t)-\Delta(t) S(t)) d t
$$

We compute the differential of the discounted portfolio value $e^{-r t} X(t)$, the differential of the discounted call price $e^{-r t} c(t, S(t))$, and set these two equal. This results in the delta-hedging rule (4.5.11),

$$
\Delta(t)=c_{x}(t, S(t)),
$$

and the Black-Scholes-Merton partial differential equation (4.5.14),

$$
c_{t}(t, x)+r x c_{x}(t, x)+\frac{1}{2} \sigma^{2} x^{2} c_{x x}(t, x)=r c(t, x) .
$$

In addition to satisfying this partial differential equation, the function $c(t, x)$ must satisfy the boundary conditions

$$
c(T, x)=(x-K)^{+}, c(t, 0)=0, \lim _{x \rightarrow \infty}\left[c(t, x)-\left(x-e^{-r(T-t)} K\right)\right]=0 .
$$

The function satisfying these conditions is (see (4.5.19))

$$
c(t, x)=x N\left(d_{+}(T-t, x)\right)-K e^{-r(T-t)} N\left(d_{-}(T-t, x)\right),
$$

where

$$
d_{\pm}(\tau, x)=\frac{1}{\sigma \sqrt{\tau}}\left[\log \frac{x}{K}+\left(r \pm \frac{\sigma^{2}}{2}\right) \tau\right] .
$$

Using the function given by (4.8.14), if one starts with initial capital $X(0)=c(0, S(0))$ and uses the delta-hedging rule $(4.8 .13)$, then at every time $t, X(t)=c(t, S(t))$. In particular, at the final time, the value of the hedging portfolio is $X(T)=c(T, S(T))=(S(T)-K)^{+}$almost surely. The short position in the European call has been hedged.

Lévy's Theorem, Theorem 4.6.4, says that if $M(t)$ is a continuous martingale starting at $M(0)=0$ and if $[M, M](t)=t$ (i.e., $d M(t) d M(t)=d t)$, then $M(t)$ is a Brownian motion. If $M_{1}(t)$ and $M_{2}(t)$ are two such processes and $\left[M_{1}, M_{2}\right](t)=0$ (i.e., $\left.d M_{1}(t) d M_{2}(t)=0\right)$, then $M_{1}(t)$ and $M_{2}(t)$ are independent Brownian motions (Theorem 4.6.5). One can use this theorem to construct independent Brownian motions from correlated Brownian motions and vice versa (see Exercise 4.13).

A Gaussian process $X(t)$ is one for which $X\left(t_{1}\right), X\left(t_{2}\right), \ldots X\left(t_{n}\right)$ are jointly normally distributed whenever $0<t_{1}<t_{2}<\cdots<t_{n}$ (Definition 4.7.1). Because the joint distribution of jointly normal random variables is determined by means, variances, and covariances, the distribution of a Gaussian process is determined by its mean function $m(t)=\mathbb{E} X(t)$ and covariance function $c(s, t)=\operatorname{Cov}(X(s), X(t))$. Brownian motion is a Gaussian process with $m(t)=0$ and $c(s, t)=s \wedge t$ (Example 4.7.2). If $\Delta(u)$ is nonrandom, then $I(t)=\int_{0}^{t} \Delta(u) d W(u)$ is a Gaussian process with $m(t)=0$ and $c(s, t)=\int_{0}^{s \wedge t} \Delta^{2}(u) d u$ (Example 4.7.3). The Brownian bridge from $a$ to $b$ on $[0, T]$ is a Gaussian process with $m(t)=\frac{(T-t) a+b t}{T}$ for $t \in[0, T]$ and $c(s, t)=s \wedge t-\frac{s t}{T}$ for $s, t \in[0, T]$ (see Subsection 10.7.2). The Brownian bridge from $a$ to $b$ on $[0, T]$ is the process one obtains by starting a Brownian motion at $a$ at time $t=0$ and conditioning on $W(T)=b$ (see Subsection 10.7.5).

\subsection{Notes}

The modern theory of stochastic calculus developed from the work of Itô [92]. Not only did Itô define the integral with respect to Brownian motion, but he also developed the change-of-variable formula commonly called Itô's rule or Itô's formula. As demonstrated in this chapter, this formula is at the heart of a wide range of useful calculations. An amazing twist to the story of stochastic calculus has recently emerged. In February 1940, the French National Academy of Sciences received a document from W. Doeblin, a French soldier on the German front. Doeblin died shortly thereafter, and the document remained sealed until May 2000 . When it was opened, the document was found to contain a construction of the stochastic integral slightly different from Itô's and a clear statement of the change-of-variable formula. Doeblin's work [52], Yor's [166] analysis of the work, and a detailed history by Bru [24] of the context of the work appeared in the December 2000 issue of Comptes Rendus de L'Académie des Sciences. An English translation of this material is [25]. Because of this remarkable development, in this text the change-ofvariable formula is called the Itô-Doeblin formula.

We have defined the Itô integral $\int_{0}^{T} \Delta^{2}(t) d W(t)$ under the condition

$$
\mathbb{E} \int_{0}^{T} \Delta^{2}(t) d t<\infty .
$$

The integral can be defined under the weaker condition

$$
\int_{0}^{T} \Delta^{2}(t)<\infty \text { almost surely }
$$

but then is not guaranteed to be a martingale. It is still a local martingale, a topic discussed in advanced books on stochastic calculus (e.g., [101]). In this text, we do not consider local martingales. We work only under the condition (4.3.1), and every Itô integral we encounter is a martingale.

Brownian motion was introduced to finance by Bachelier [6]. Samuelson [143], [145] presents the argument that geometric Brownian motion is a good model for stock prices. The application of stochastic calculus to finance began with the work of Merton [121]. (The paper [121] and many other papers by Merton that use stochastic calculus in finance are collected in Merton [124].) The Black-Scholes-Merton formula is based on the geometric Brownian motion model for stock prices. However, no-arbitrage pricing theory has now moved far beyond this assumption. As seen in this and subsequent chapters, this theory and the accompanying risk-neutral pricing formula can be applied in the presence of a time-varying random volatility, a time-varying random mean rate of return, and a time-varying random interest rate.

Many finance books, including (in order of increasing mathematical difficulty) Hull [87], Dothan [54], and Duffie [56], include sections on Itô's integral and the Itô-Doeblin formula. Some other books on dynamic models in finance are Cox and Rubinstein [43], Huang and Litzenberger [86], Ingersoll [91], and Jarrow [97]. A comprehensive text is Wilmott [164]. Some good references for practitioners are Baxter and Rennie [8] (reviewed in [134]), Björk [11] (reviewed in [135]), and Musiela and Rutkowski [126] (reviewed in [134]). More mathematical texts on stochastic calculus with applications to finance are Lamberton and Lapeyre [105] (reviewed in [134]) and Steele [150] (reviewed in [136]). Other texts on stochastic calculus are Chung and Williams [36], Karatzas and Shreve [101], Øksendal [129], and Protter [133]. Karatzas and Shreve [102] is a sequel to [101] that focuses on finance. Protter [133] is the easiest place to learn about stochastic calculus for processes with jumps, and this is not at all easy. We introduce this topic in Chapter 11.

No-arbitrage pricing theory and the accompanying risk-neutral pricing formula is predicated on the assumption that there is no arbitrage in the market. An arbitrage is defined to be a trading strategy which begins with zero capital and at a later time has positive capital with positive probability without having any risk of loss. Absence of arbitrage is similar to but different from the efficient market hypothesis, which asserts that technical analysis of stock prices is of no value. This hypothesis asserts that patterns in stock prices may be useful to estimate the parameters of the distribution of future returns, but they do not provide clues to whether the next price movement will be up or down. In particular, technical analysis does not permit one to outperform the market. This hypothesis could be violated in a way which permits one to outperform the market with high probability without actually admitting arbitrage because there is still a nonzero probability of underperforming the market. This is sometimes called statistical arbitrage. An empirical study supporting the efficient market hypothesis is Fama [64], which also discusses distributions that fit stock prices better than geometric Brownian motion. A criticism of the efficient market hypothesis is provided by LeRoy [106], and a recent paper that finds long-range dependence (but not much) in stock price data is Willinger, Taqqu, and Teverovsky [163]. A provocative article on the source of stock price movements is Black [14].

Geometric fractional Brownian motion has been proposed as an alternative model for stock prices because it has fatter tails than geometric Brownian motion. One can assume such a model and compute the prices of derivative 
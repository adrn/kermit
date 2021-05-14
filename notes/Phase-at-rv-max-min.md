# Reparametrizing the Two Body Problem

When fitting radial velocity (RV) data for exoplanet or binary star systems, the
RV model is typically expressed as:
$$
    v_r(t) = v_0 + K\,\left[\cos(\omega + f(t)) + e\,\cos(\omega)\right]
$$
where $v_0$ is the systemic (barycentric) velocity of the system relative to the
Sun, $K$ is the velocity semi-amplitude, $\omega$ is the argument of pericenter,
and $f(t)$ is the True Anomaly. The True Anomaly is computed from the Mean
Anomaly $M(t)$
$$
    M = \frac{2\pi}{P} \, (t - t_{\mathrm{peri}})
$$
through a transcendental equation that relates the Mean Anomaly to the Eccentric
Anomaly $E(t)$:
$$
    M = E - e\,\sin E
$$
The Eccentric Anomaly is related to the True Anomaly by:
$$
\begin{aligned}
    \tan\frac{f}{2} &= \sqrt{\frac{1+e}{1-e}} \, \tan\frac{E}{2} \\
    \tan\frac{E}{2} &= \sqrt{\frac{1-e}{1+e}} \, \tan\frac{f}{2}
\end{aligned}
$$

In the above parametrization, the full set of parameters is $(v_0, K, P, e,
\omega, t_{\mathrm{peri}})$.

Another common representation of the Mean Anomaly is in terms of the phase of
pericenter $M_0$, which specifies the phase at which pericenter occurs in the
angle defined by time relative to a reference time $t_{\mathrm{ref}}$, which is
often taken to be the minimum or mean observation time. So, an alternate
definition of Mean Anomaly is:
$$
    M = \frac{2\pi}{P} \, (t - t_{\mathrm{ref}}) - M_0
$$
With this convention, $t_{\mathrm{peri}}$ is replaced by $M_0$ as a parameter so
the list of parameters is $(v_0, K, P, e, \omega, M_0)$, but it is
straightforward to transform between $M_0$ and $t_{\mathrm{peri}}$. In either
case, the Mean Anomaly is defined such at $M=0$ (mod $2\pi$) occurs at
pericenter.

TODO: discussion about how this is not a good parametrization for MCMC sampling,
and history (Ford 2005).

## Replacing $\omega, M_0$ with phase of max/min velocity $M^*_{\mathrm{max}}, M^*_{\mathrm{min}}$

$$
\begin{aligned}
M^*_{\mathrm{max}} &= \frac{2\pi}{P} \, (t_{\mathrm{max}} - t_{\mathrm{ref}}) \\
M^*_{\mathrm{min}} &= \frac{2\pi}{P} \, (t_{\mathrm{min}} - t_{\mathrm{ref}})
\end{aligned}
$$
Note: $M^*_{\mathrm{max}}, M^*_{\mathrm{min}}$ are not actually mean anomalies!
They are mean anomaly $+ M_0$, because what we observe are phases relative to
the reference time $t_{\mathrm{ref}}$.

So, how do we relate these quantities to $(\omega, M_0)$?

As far as I can figure, there is no direct, closed-form transformation, but I
have figured out a procedure to do the transformation in terms of a different
transcendental equation. This method is built on the intuition (guess) that the
quantity
$$
    \Delta M = M^*_{\mathrm{min}} - M^*_{\mathrm{max}}
$$
will be related to the argument of pericenter. Defining also
$$
\begin{aligned}
    \Delta E &= E_{\mathrm{min}} - E_{\mathrm{max}} \\
    \Sigma E &= E_{\mathrm{min}} + E_{\mathrm{max}}
\end{aligned}
$$
we can re-express the classic Kepler transcendental equation as
$$
\begin{aligned}
    \Delta M &= \Delta E - e\, (\sin E_{\mathrm{min}} - \sin E_{\mathrm{max}})\\
    \frac{\Delta M}{2} &= \frac{\Delta E}{2} - e\, \sin\frac{\Delta E}{2} \, \cos\frac{\Sigma E}{2}
\end{aligned}
$$
where the last line makes use of the "difference of sines" trig identity.

Using the relationship between $E$ and $f$,
$$
    (E_{\mathrm{min}} - E_{\mathrm{max}})/2 =
        \arctan\left[ \sqrt{\frac{1-e}{1+e}} \,
            \tan \frac{f_{\mathrm{min}}}{2} \right] -
        \arctan\left[ \sqrt{\frac{1-e}{1+e}} \,
            \tan \frac{f_{\mathrm{max}}}{2} \right]
$$
However, we know the values of $f_{\mathrm{min}}$ and $f_{\mathrm{max}}$: From
Equation (1), we know that the RV is maximum when $f=f_{\mathrm{max}} \equiv
-\omega$, and minimum when $f=f_{\mathrm{min}} \equiv \pi-\omega$.
Plugging in these values and doing some trig algebra (using the difference of
arctan's trig identity), this simplifies to
$$
    \Delta E / 2 = \arctan\left(\frac{\sqrt{1 - e^2}}{e\,\sin\omega}\right)
$$
and similarly, for the sum of $E_{\mathrm{min}} + E_{\mathrm{max}}$,
$$
    \Sigma E / 2 = \arctan\left(\frac{\sqrt{1 - e^2}}{\tan\omega}\right)
        \quad .
$$
So, the relationship between $\omega$, $e$, and $\Delta M$ is defined by the
expression
$$
\begin{aligned}
    \frac{\Delta M}{2} &= \frac{\Delta E}{2} -
        e\, \sin\frac{\Delta E}{2} \, \cos\frac{\Sigma E}{2} \\
    \frac{\Delta M}{2} &=
        \arctan\left(\frac{\sqrt{1 - e^2}}{e\,\sin\omega}\right) -
        \,e \,
            \sin\left(\arctan\left(\frac{\sqrt{1 - e^2}}{e\,\sin\omega}\right)\right) \,
            \cos\left(\arctan\left(\frac{\sqrt{1 - e^2}}{\tan\omega}\right)\right)
\end{aligned}
$$
and the phase $M_0$ can be retrieved with
$$
\begin{aligned}
    M_0 &= M^*_{\mathrm{max}} - M_{\mathrm{max}} \\
    M_{\mathrm{max}} &= E_{\mathrm{max}} - e \, \sin E_{\mathrm{max}} \\
    E_{\mathrm{max}} &= 2\,\arctan\left(\sqrt{\frac{1-e}{1+e}} \,
        \tan\frac{f_{\mathrm{max}}}{2}\right) \\
    &= -2\,\arctan\left(\sqrt{\frac{1-e}{1+e}} \,
        \tan\frac{\omega}{2}\right) \\
\end{aligned}
$$



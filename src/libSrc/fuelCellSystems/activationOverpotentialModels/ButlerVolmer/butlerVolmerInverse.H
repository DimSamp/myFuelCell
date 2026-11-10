/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by the original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

InNamespace
    Foam::activationOverpotentialModels

Description
    Inverse of the dimensionless Butler-Volmer relation

    \f[
        h(y) = \exp(\alpha y) - \exp(-(1 - \alpha) y) = r ,
        \qquad 0 < \alpha < 1
    \f]

    with y = n F eta/(R T) and r = (current density)/(kinetic prefactor),
    both positive in the forward direction of the electrode reaction.

    h is strictly increasing, so the root is unique. It is found by
    Newton-Raphson, safeguarded by a bracket: the root lies between 0 and
    the Tafel limit of the dominant branch,

        r > 0:   0 < y <= ln(1 + r)/alpha
        r < 0:   -ln(1 - r)/(1 - alpha) <= y < 0

    and a Newton step that leaves the current bracket is replaced by a
    bisection step. The iteration stops when the last step is below
    tol*|y|; with the quadratic convergence of Newton-Raphson the root is
    then accurate to round-off.

\*---------------------------------------------------------------------------*/

#ifndef butlerVolmerInverse_H
#define butlerVolmerInverse_H

#include "scalar.H"
#include "label.H"
#include <cmath>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace activationOverpotentialModels
{

//- Solve h(y) = r. Returns y; sets dhdy = h'(y) > 0 and the number of
//- iterations nIter (nIter > maxIter: not converged)
inline scalar butlerVolmerInverse
(
    const scalar r,
    const scalar alpha,
    scalar& dhdy,
    label& nIter,
    const scalar tol = 1e-12,
    const label maxIter = 100
)
{
    const scalar beta = scalar(1) - alpha;

    //- A non-finite current is passed on, flagged as not converged
    if (!std::isfinite(r))
    {
        dhdy = r;
        nIter = maxIter + 1;
        return r;
    }

    //- Bracket of the root and start value: the small-current limit y = r
    //- (h'(0) = 1) or the Tafel limit, whichever is closer to zero
    scalar yLo = 0;
    scalar yHi = 0;
    scalar y = 0;

    if (r > 0)
    {
        yHi = std::log1p(r)/alpha;
        y = (r < yHi ? r : yHi);
    }
    else if (r < 0)
    {
        yLo = -std::log1p(-r)/beta;
        y = (r > yLo ? r : yLo);
    }

    for (nIter = 1; nIter <= maxIter; ++nIter)
    {
        //- expm1 keeps the difference accurate for |y| << 1
        const scalar ea = std::expm1(alpha*y);
        const scalar eb = std::expm1(-beta*y);

        const scalar res = (ea - eb) - r;
        dhdy = alpha*(ea + 1) + beta*(eb + 1);

        //- h is increasing: tighten the bracket
        if (res > 0)
        {
            yHi = y;
        }
        else
        {
            yLo = y;
        }

        scalar yNew = y - res/dhdy;

        if (!(yNew >= yLo && yNew <= yHi))
        {
            yNew = 0.5*(yLo + yHi);
        }

        const scalar dy = yNew - y;
        y = yNew;

        if (std::abs(dy) <= tol*std::abs(y))
        {
            break;
        }
    }

    dhdy = alpha*std::exp(alpha*y) + beta*std::exp(-beta*y);

    return y;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace activationOverpotentialModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

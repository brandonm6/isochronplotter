"""
Copyright (c) 2018, Lawrence Livermore National Security, LLC. Produced at the Lawrence Livermore National Laboratory.
Written by Reto Trappitsch, trappitsch1@llnl.gov

LLNL-CODE-745740 All rights reserved. This file is part of MahonFitting v1.0

Please also read this link - Our Disclaimer (https://github.com/LLNL/MahonFitting/blob/master/DISCLAIMER) and
GNU General Public License (https://github.com/LLNL/MahonFitting/blob/master/LICENSE).

This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public
License (as published by the Free Software Foundation) version 2, dated June 1991.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the
IMPLIED WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the terms and conditions of the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program; if not, write to the
Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
"""
# calculate linear regression w/ error calculation according Mahon 1996 (New York regression)


# Streamlined to remove GUI and require dataframe as input
import numpy as np


class Mahon:
    def __init__(self, df, sigma_unc=1):
        """
        Calculates linear regression w/ error calculation according Mahon 1996 (New York regression)
        :param df: dataframe with the following columns:
                   [x], [x unc], [y], [y unc], [corr]
                    - [x] are the x values
                    - [x unc] are the uncertainties for x
                    - [y] are the y values
                    - [y unc] are the uncertainties for y
                    - [corr] are the correlation values

        :param sigma_unc: number of sigma uncertainties (either 1 or 2)
        """

        # some variables you can call
        self.slope = None
        self.yinter = None
        self.xinter = None
        self.slopeunc = None
        self.yinterunc = None
        self.xinterunc = None
        self.mswd = None

        # some variables to remember
        self.xbar = None
        self.ybar = None

        # how many sigma uncertainties are in the file?
        self.sigvar = sigma_unc

        # create variables for holding the data
        self.xdat = df["x"].to_numpy()
        self.ydat = df["y"].to_numpy()
        self.xunc = (df["x unc"] / self.sigvar).to_numpy()
        self.yunc = (df["y unc"] / self.sigvar).to_numpy()
        self.p = df["corr"].to_numpy()

        # run the calculation with the loaded x and y data
        self.calcparams()
        self.calcunc()
        self.calcunc(calcxintunc=True)

        # calculate the MSWD
        self.calcmswd()

    def calcparams(self):
        """
        Calculate the parameters, both intercepts as well as the slope of the regression
        :return:
        """
        bcalc = 1  # initial guess for slope
        bold = 0  # to compare later

        # read in selfs
        xdat = self.xdat
        xunc = self.xunc
        ydat = self.ydat
        yunc = self.yunc

        # run thorough the while loop
        whilecounter = 0
        whilecountermax = 1e5
        while np.abs((bold - bcalc) / bcalc) > 1.e-10 and whilecounter < whilecountermax:
            whilecounter += 1
            # prep for while loop, start before this line and compare bold to bcalc
            bold = bcalc
            # calculate xbar and ybar
            xbar = 0
            ybar = 0
            weightsum = 0
            for it in range(len(xdat)):
                wi = calc_wi(xunc[it], yunc[it], bcalc, self.p[it])
                xbar += xdat[it] * wi
                ybar += ydat[it] * wi
                weightsum += wi
            xbar /= weightsum
            ybar /= weightsum

            # now calculate b
            btop = 0  # top sum
            bbot = 0  # bottom sum

            for it in range(len(xdat)):
                xi = xdat[it]
                yi = ydat[it]
                sxi = xunc[it]
                syi = yunc[it]
                pit = self.p[it]
                wi = calc_wi(sxi, syi, bcalc, pit)
                ui = xi - xbar
                vi = yi - ybar
                # add to sums
                btop += wi ** 2. * vi * (ui * syi ** 2. + bcalc * vi * sxi ** 2. - pit * vi * sxi * syi)
                bbot += wi ** 2. * ui * (ui * syi ** 2. + bcalc * vi * sxi ** 2. - bcalc * pit * ui * sxi * syi)

            # new bcalc
            bcalc = btop / bbot

        # error message if whilecounter timed out
        if whilecounter == whilecountermax:
            print('Warning! Your calculation might not have converged ' +
                  'properly. The difference between the last calculated '
                  'slope and the current slope is: ' + str(np.abs((bold - bcalc) / bcalc)) +
                  ' You can ignore this ' + 'message if this is an acceptable convergence for you.')

        # now that slope is determined, calculate the y intercept
        self.yinter = ybar - bcalc * xbar

        # now done, so write back slope
        self.slope = bcalc

        # calculate x intercept
        self.xinter = -self.yinter / self.slope

        # write back xbar and ybar
        self.xbar = xbar
        self.ybar = ybar

    def calcunc(self, calcxintunc=False):
        """
        Calculates the uncertainty of the slope and y
        :param calcxintunc: If it needs to calculate the x uncertainty, then set this to true
        :return:
        """
        if calcxintunc:
            # read in selfs
            # since this is for x uncertainty, simply switch it x and y.
            xdat = self.ydat
            xunc = self.yunc
            ydat = self.xdat
            yunc = self.xunc
            xbar = self.ybar
            ybar = self.xbar
            b = 1. / self.slope
        else:
            # read in selfs
            xdat = self.xdat
            xunc = self.xunc
            ydat = self.ydat
            yunc = self.yunc
            xbar = self.xbar
            ybar = self.ybar
            b = self.slope

        # let us first calculate the derivatives
        # dell theta / dell b (dthdb) calculation
        sum1 = 0.
        sum2 = 0.
        for it in range(len(xdat)):
            xi = xdat[it]
            yi = ydat[it]
            sxi = xunc[it]
            syi = yunc[it]
            pit = self.p[it]
            wi = calc_wi(xunc[it], yunc[it], b, pit)
            ui = xi - xbar
            vi = yi - ybar
            sxyi = pit * sxi * syi
            sum1 += wi ** 2. * (
                    2 * b * (ui * vi * sxi ** 2. - ui ** 2. * sxyi) + (ui ** 2. * syi ** 2. - vi ** 2 * sxi ** 2.))
            sum2 += wi ** 3. * (sxyi - b * sxi ** 2.) * (b ** 2. * (ui * vi * sxi ** 2 - ui ** 2 * sxyi) +
                                                         b * (ui ** 2 * syi ** 2 - vi ** 2 * sxi ** 2) -
                                                         (ui * vi * syi ** 2 - vi ** 2 * sxyi))
        dthdb = sum1 + 4. * sum2

        # calculate the sum of all weights
        wksum = 0.
        for it in range(len(xdat)):
            wksum += calc_wi(xunc[it], yunc[it], b, self.p[it])

        # now calculate sigasq and sigbsq
        sigasq = 0.
        sigbsq = 0.
        for it in range(len(xdat)):
            sxi = xunc[it]
            syi = yunc[it]
            pit = self.p[it]
            wi = calc_wi(sxi, syi, b, pit)
            sxyi = pit * sxi * syi

            # calculate dell theta / dell xi and dell theta / dell yi
            dthdxi = 0.
            dthdyi = 0.
            for jt in range(len(xdat)):
                xj = xdat[jt]
                yj = ydat[jt]
                sxj = xunc[jt]
                syj = yunc[jt]
                pjt = self.p[jt]
                wj = calc_wi(sxj, syj, b, pjt)
                uj = xj - xbar
                vj = yj - ybar
                sxyj = pjt * sxj * syj
                # add to dthdxi and dthdyi
                dthdxi += wj ** 2. * (kron(it, jt) - wi / wksum) * (b ** 2 * (vj * sxj ** 2 - 2 * uj * sxyj) +
                                                                    2 * b * uj * syj ** 2 - vj * syj ** 2)
                # correct equation! not equal to equation 21 in Mahon (1996)
                dthdyi += wj ** 2. * (kron(it, jt) - wi / wksum) * (b ** 2 * uj * sxj ** 2 + 2 * vj * sxyj -
                                                                    2 * b * vj * sxj ** 2. - uj * syj ** 2)

            # now calculate dell a / dell xi and dell a / dell yi
            dadxi = -b * wi / wksum - xbar * dthdxi / dthdb
            dadyi = wi / wksum - xbar * dthdyi / dthdb

            # now finally add to sigasq and sigbsq
            sigbsq += dthdxi ** 2. * sxi ** 2. + dthdyi ** 2. * syi ** 2. + 2 * sxyi * dthdxi * dthdyi
            sigasq += dadxi ** 2. * sxi ** 2. + dadyi ** 2. * syi ** 2. + 2 * sxyi * dadxi * dadyi

        # now divide sigbsq
        sigbsq /= dthdb ** 2.

        # now write slope uncertainty and y intercept uncertainty back to class
        if calcxintunc:
            self.xinterunc = np.sqrt(sigasq) * float(self.sigvar)
        else:
            self.yinterunc = np.sqrt(sigasq) * float(self.sigvar)
            self.slopeunc = np.sqrt(sigbsq) * float(self.sigvar)

    def calcmswd(self):
        xdat, ydat, xunc, yunc = self.xdat, self.ydat, self.xunc, self.yunc
        mswd = 0.
        for it in range(len(xdat)):
            xi = xdat[it]
            yi = ydat[it]
            sxi = xunc[it]
            syi = yunc[it]
            pit = self.p[it]
            wi = calc_wi(sxi, syi, self.slope, pit)
            mswd += wi * (yi - self.slope * xi - self.yinter) ** 2.

        # now divide by degrees of freedom minus 2, since 2 fixed parameters
        mswd /= (len(xdat) - 2.)
        self.mswd = mswd


def calc_wi(sx, sy, b, p):
    return 1. / (sy ** 2 + b ** 2 * sx ** 2 - 2 * b * p * sx * sy)


def kron(i, j):
    # calculates Kronecker delta
    if i == j:
        return 1.
    else:
        return 0.

*! oaxaca_twofold_lee_v2 v2.9  -- 2026-08-23
*! ============================================================================
*! TWOFOLD OAXACA-BLINDER DECOMPOSITION WITH A LEE (1983) POLYCHOTOMOUS
*! SELECTION CORRECTION, REPORTED UNDER THE NONCOMMITTAL SCHEME OF NEUMAN AND
*! OAXACA (2004).
*!
*! Authors: Mauricio Quinones-Dominguez, Sergio A. Barona-Montoya
*! Copyright (c) 2026. Released under the MIT License; see LICENSE.
*!
*! ----------------------------------------------------------------------------
*! WHAT IT REPORTS
*! ----------------------------------------------------------------------------
*!   Difference   observed gap, Ybar(ref) - Ybar(other)
*!   Endowments   (Xbar_ref - Xbar_other)' beta_ref
*!   Returns      Xbar_other' (beta_ref - beta_other)
*!   Selectivity  theta_ref*hbar_ref - theta_other*hbar_other
*!   WageOffers   Xbar_ref'beta_ref - Xbar_other'beta_other
*!
*!   The first four components are the noncommittal decomposition of Neuman
*!   and Oaxaca (2004) and add up exactly to the observed gap. WageOffers is
*!   the estimand obtained by netting the selectivity term out of the observed
*!   gap: the gap that would prevail in the absence of differential selection.
*!
*! ----------------------------------------------------------------------------
*! VERSION HISTORY
*! ----------------------------------------------------------------------------
*! v2.9  FIX: Stata prefixes the names of coefficients dropped for collinearity
*!       with "o.", so e(b) may contain "o.sector_5" instead of "sector_5".
*!       Matching on the exact name failed with rc=111. Because the set of
*!       collinear regressors differs across cells, the failure appeared only
*!       in some of them. Any factor prefix ("o.", "b.", "1b.", "2o.") is now
*!       stripped before the comparison. A dropped coefficient equals zero and
*!       contributes zero to both components, which is the correct treatment;
*!       including it only keeps the dimensions of beta and Xbar aligned.
*! v2.8  Exact variance for WageOffers and a second closure check, e(closure2).
*! v2.7  Fixes the r(111) error caused by reading tempnames after -ereturn post-.
*! v2.6  Standard error for Difference; vrobust option; diagnostics.
*! v2.5  Estimation sample rebuilt with markout.
*! v2.4  Option nowls: exact accounting closure.
*! v2.3  Correction terms located through the _m prefix.
*! v2.2  The first stage retains all three states.
*!
*! ----------------------------------------------------------------------------
*! SYNTAX
*! ----------------------------------------------------------------------------
*!   oaxaca_twofold_lee_v2 depvar indepvars [if] [in] [aw pw fw], ///
*!       by(group) select(zvarlist) [ref(#) status(name) genpfx(name)
*!        innerboot(#) nowls vrobust nodetail]
*!
*!   depvar   MUST be missing outside the sector of interest. Do not pass
*!            -if status_var==#-: the first stage needs all three states.
*!   nowls    second stage by OLS. RECOMMENDED: it makes the closure exact.
*!   vrobust  heteroskedasticity-robust variance; implies nowls.
*!
*!   Requires selmlog (Bourguignon, Fournier and Gurgand).
*!
*!   Detailed components are not comparable across cells when the set of
*!   regressors dropped for collinearity differs between them. The aggregate
*!   components are comparable, because the explained term sums over all
*!   dummies and is invariant to which one is omitted.
*! ============================================================================

program define oaxaca_twofold_lee_v2, eclass
    version 16

    // ======================================================================
    // 0. SYNTAX AND OPTIONS
    // ======================================================================

    syntax varlist(min=2) [if] [in] [aweight pweight fweight],           ///
        BY(varname) SELect(varlist)                                      ///
        [ REF(string) STATus(name) GENpfx(name)                          ///
          INNERboot(integer 0) NOWLS VROBust NODEtail ]

    gettoken depvar indepvars : varlist
    local x_ols "`indepvars'"

    if "`status'" == "" local status status_var
    if "`genpfx'" == "" local genpfx _h

    // vrobust takes b and V from an OLS refit, so the second stage must be
    // OLS for the reported coefficients to be those of selmlog.
    if "`vrobust'" != "" & "`nowls'" == "" {
        di as txt "Note: vrobust implies an OLS second stage; nowls assumed."
        local nowls "nowls"
    }
    local wlsopt = cond("`nowls'" == "", "wls", "")

    // touse must NOT drop observations because `depvar' is missing: those
    // observations identify the first stage. The sector restriction travels
    // through the dependent variable, not through the sample.
    marksample touse, novarlist
    markout `touse' `by' `select' `status'

    // ======================================================================
    // 1. GROUPS AND REFERENCE STRUCTURE
    // ======================================================================

    qui levelsof `by' if `touse', local(glevels)
    local ng : word count `glevels'
    if `ng' != 2 {
        di as err "`by' must take exactly two values in the sample."
        exit 198
    }
    local lev1 : word 1 of `glevels'
    local lev2 : word 2 of `glevels'

    if "`ref'" == "" {
        local gm `lev2'
        di as txt "Note: ref() not specified; using `by' == `gm'."
    }
    else {
        local gm `ref'
        if (`gm' != `lev1') & (`gm' != `lev2') {
            di as err "ref(`gm') is not a value of `by'."
            exit 198
        }
    }
    if `gm' == `lev1'  local gf `lev2'
    else               local gf `lev1'

    local lblname : value label `by'
    local gm_lbl "`gm'"
    local gf_lbl "`gf'"
    if "`lblname'" != "" {
        local gm_lbl : label `lblname' `gm'
        local gf_lbl : label `lblname' `gf'
    }
    di as txt _n "Reference structure : `by' == `gm' [`gm_lbl']"
    di as txt    "Comparison group    : `by' == `gf' [`gf_lbl']"
    di as txt    "Reported gap        : Ybar(`gm_lbl') - Ybar(`gf_lbl')"
    di as txt    "Second stage        : " ///
                 cond("`nowls'" == "", "WLS (selmlog default)", "OLS")
    di as txt    "Variance            : " ///
                 cond("`vrobust'" == "", "e(V) from selmlog, not robust", ///
                                         "heteroskedasticity-robust") _n

    local bootopt ""
    if `innerboot' > 0 local bootopt "bootstrap(`innerboot')"

    tempvar one
    qui gen byte `one' = 1

    // ======================================================================
    // 2. ESTIMATION BY GROUP
    // ======================================================================
    // Suffix "m" = reference group, "f" = comparison group.
    // ======================================================================

    foreach s in m f {

        if "`s'" == "m"  local gval `gm'
        else             local gval `gf'

        //------------------------------------------------------------------
        // 2.1 Check: the first stage needs all three states
        //------------------------------------------------------------------
        qui levelsof `status' if `touse' & `by' == `gval', local(_stlev)
        local nst : word count `_stlev'
        if `nst' < 3 {
            di as err _n "`status' takes only `nst' value(s) in this sample."
            di as err "The Lee correction needs all three states in the FIRST"
            di as err "stage. Do not pass -if `status'==#-. Restrict the"
            di as err "sector through the dependent variable instead:"
            di as err "    gen dep_formal = dep_var if `status'==1"
            exit 148
        }

        qui count if `touse' & `by' == `gval' & !missing(`depvar')
        if r(N) < 30 {
            di as err "Only `r(N)' observations with `depvar' nonmissing" ///
                      " for `by' == `gval'."
            exit 2001
        }

        //------------------------------------------------------------------
        // 2.2 Two-step estimation (selmlog)
        //------------------------------------------------------------------
        // The selmlog warning about observations with missing values in the
        // main equation reports the intended behaviour: it is what makes the
        // correction identifiable.
        //------------------------------------------------------------------
        capture drop `genpfx'*
        di as txt "Estimating for `by' == `gval' ..."

        selmlog `depvar' `x_ols' if `touse' & `by' == `gval',             ///
            select(`status' = `select') lee `wlsopt' `bootopt'            ///
            mloptions(level(95)) gen(`genpfx')

        tempname bfull Vfull
        matrix `bfull' = e(b)
        matrix `Vfull' = e(V)
        local cn : colnames `bfull'
        local ce : coleq    `bfull'
        local k = colsof(`bfull')

        //------------------------------------------------------------------
        // 2.3 Locating the correction terms
        //------------------------------------------------------------------
        // selmlog names the COEFFICIENT _m<j> inside the wage equation and
        // places Sigma<j> and rho<j> in a separate equation called
        // "Anciliary". The corresponding VARIABLE in the dataset is named
        // `genpfx'<j>, through the gen() option. They are linked by the
        // numeric suffix: the coefficient gives theta_j, the variable hbar_j.
        //------------------------------------------------------------------
        local hvars ""
        local rho_`s' = .
        forvalues c = 1/`k' {
            local nm : word `c' of `cn'
            local eq : word `c' of `ce'

            if inlist(lower("`eq'"), "anciliary", "ancillary") {
                if substr(lower("`nm'"), 1, 3) == "rho" ///
                    local rho_`s' = `bfull'[1, `c']
                continue
            }
            local nm_lim = regexr("`nm'", "^[0-9]*[bo]?\.", "")
            if substr("`nm_lim'", 1, 2) == "_m" {
                local jsuf = substr("`nm_lim'", 3, .)
                capture confirm numeric variable `genpfx'`jsuf'
                if _rc {
                    di as err "Found coefficient `nm' but variable" ///
                              " `genpfx'`jsuf' does not exist."
                    di as err "Check the prefix passed in genpfx()."
                    exit 111
                }
                local hvars "`hvars' `genpfx'`jsuf'"
            }
        }
        if "`hvars'" == "" {
            di as err _n "No correction terms (_m*) found in e(b)."
            di as err "Column names: `cn'"
            exit 111
        }
        di as txt "  correction term(s): `hvars'"

        //------------------------------------------------------------------
        // 2.4 Wage-equation sample
        //------------------------------------------------------------------
        // e(sample) from selmlog cannot be used: under -nowls- it comes back
        // empty. It is rebuilt by reproducing the listwise deletion of the
        // second step. The count printed below must match the "Number of obs"
        // of the second-stage regression.
        //------------------------------------------------------------------
        tempvar es`s'
        qui gen byte `es`s'' = (`touse' & `by' == `gval' & !missing(`depvar'))
        markout `es`s'' `x_ols' `hvars'
        qui count if `es`s''
        local nobs_`s' = r(N)
        di as txt "  wage-equation observations: `nobs_`s''"
        if `nobs_`s'' == 0 {
            di as err "Empty sample for `by' == `gval'."
            exit 2000
        }

        //------------------------------------------------------------------
        // 2.5 Robust refit (optional)
        //------------------------------------------------------------------
        // Under nowls the second stage is an OLS regression of depvar on
        // x_ols and the correction terms. Refitting with -regress, robust-
        // reproduces the selmlog coefficients exactly and only replaces the
        // variance matrix. This matters because the uncorrected branch of the
        // analysis uses robust standard errors: without it, the two panels of
        // each figure are not comparable in terms of inference.
        //------------------------------------------------------------------
        if "`vrobust'" != "" {
            qui regress `depvar' `x_ols' `hvars' if `es`s'', robust
            matrix `bfull' = e(b)
            matrix `Vfull' = e(V)
            local cn : colnames `bfull'
            local ce : coleq    `bfull'
            local k = colsof(`bfull')
            local hnames "`hvars'"      // regress names them by variable
        }
        else {
            local hnames ""             // selmlog names them _m<j>
        }

        //------------------------------------------------------------------
        // 2.6 Positions of the covariates
        //------------------------------------------------------------------
        // NOTE ON THE "o." PREFIX: Stata prepends it to the names of
        // coefficients dropped for collinearity, so e(b) may contain
        // "o.sector_5" instead of "sector_5". Matching on the exact name
        // failed with rc=111, and because the set of collinear regressors
        // differs across cells, the failure appeared only in some of them.
        //
        // regexr strips any factor prefix: "o.", "b.", "1b.", "2o.". A dropped
        // coefficient equals zero and contributes zero to both endowments and
        // returns, which is the correct treatment; including it only keeps the
        // dimensions of beta and Xbar aligned. The accounting closure
        // confirms this.
        //------------------------------------------------------------------
        local xpos ""
        local nomit 0
        foreach v in `x_ols' _cons {
            local hit = 0
            forvalues c = 1/`k' {
                local nm : word `c' of `cn'
                local eq : word `c' of `ce'
                if inlist(lower("`eq'"), "anciliary", "ancillary") continue
                local nm_lim = regexr("`nm'", "^[0-9]*[bo]?\.", "")
                if ("`nm_lim'" == "`v'") & (`hit' == 0) {
                    local xpos "`xpos'`c' "
                    local hit = 1
                    if "`nm'" != "`nm_lim'" local ++nomit
                }
            }
            if `hit' == 0 {
                di as err "Coefficient `v' not found."
                di as err "Columns: `cn'"
                exit 111
            }
        }
        if `nomit' > 0 {
            di as txt "  covariates dropped for collinearity: `nomit'" ///
                      " (zero coefficient, zero contribution)"
        }

        //------------------------------------------------------------------
        // 2.7 Positions of the correction terms
        //------------------------------------------------------------------
        local tpos ""
        if "`hnames'" != "" {
            // vrobust route: search by variable name
            foreach v of local hnames {
                local hit = 0
                forvalues c = 1/`k' {
                    local nm : word `c' of `cn'
                    local nm_lim = regexr("`nm'", "^[0-9]*[bo]?\.", "")
                    if ("`nm_lim'" == "`v'") & (`hit' == 0) {
                        local tpos "`tpos'`c' "
                        local hit = 1
                    }
                }
                if `hit' == 0 {
                    di as err "Term `v' not found after the refit."
                    di as err "Columns: `cn'"
                    exit 111
                }
            }
        }
        else {
            // non-vrobust route: search by the _m prefix of the coefficient
            forvalues c = 1/`k' {
                local nm : word `c' of `cn'
                local eq : word `c' of `ce'
                if inlist(lower("`eq'"), "anciliary", "ancillary") continue
                local nm_lim = regexr("`nm'", "^[0-9]*[bo]?\.", "")
                if substr("`nm_lim'", 1, 2) == "_m" local tpos "`tpos'`c' "
            }
        }

        //------------------------------------------------------------------
        // 2.8 Extracting the submatrices
        //------------------------------------------------------------------
        // Two blocks are taken from e(b) and e(V): covariates (beta, V_beta)
        // and correction terms (theta, V_theta). They are built element by
        // element because Stata cannot index a matrix by a vector of
        // positions.
        //
        // The covariance between beta and theta (the off-diagonal blocks) is
        // discarded. This is harmless: no individual component depends on
        // both.
        //------------------------------------------------------------------
        tempname bsub Vsub tsub Vtsub

        // --- beta and its variance ---
        local nx : word count `x_ols' _cons
        matrix `bsub' = J(1, `nx', .)
        local j = 0
        foreach c of local xpos {
            local ++j
            matrix `bsub'[1, `j'] = `bfull'[1, `c']
        }
        matrix __b_`s' = `bsub'

        matrix `Vsub' = J(`nx', `nx', 0)
        local r = 0
        foreach cr of local xpos {
            local ++r
            local q = 0
            foreach cc of local xpos {
                local ++q
                matrix `Vsub'[`r', `q'] = `Vfull'[`cr', `cc']
            }
        }
        matrix __V_b`s' = `Vsub'

        // --- theta and its variance ---
        local nt : word count `hvars'
        matrix `tsub' = J(1, `nt', .)
        local j = 0
        foreach c of local tpos {
            local ++j
            matrix `tsub'[1, `j'] = `bfull'[1, `c']
        }
        matrix __th_`s' = `tsub'

        matrix `Vtsub' = J(`nt', `nt', 0)
        local r = 0
        foreach cr of local tpos {
            local ++r
            local q = 0
            foreach cc of local tpos {
                local ++q
                matrix `Vtsub'[`r', `q'] = `Vfull'[`cr', `cc']
            }
        }
        matrix __V_th`s' = `Vtsub'

        //------------------------------------------------------------------
        // 2.9 Means over the wage-equation sample
        //------------------------------------------------------------------
        // -mean- returns the vector of means in e(b) and its
        // variance-covariance matrix in e(V), which is exactly V(Xbar),
        // respecting the weights. The mean of the dependent variable supplies
        // Ybar and its variance, which is what the Difference component
        // needs.
        //------------------------------------------------------------------
        qui mean `x_ols' `one' if `es`s'' [`weight'`exp']
        matrix __xbar_`s'  = e(b)
        matrix __Vxbar_`s' = e(V)

        qui mean `hvars' if `es`s'' [`weight'`exp']
        matrix __hbar_`s'  = e(b)
        matrix __Vhbar_`s' = e(V)

        qui mean `depvar' if `es`s'' [`weight'`exp']
        matrix __ybar_`s'  = e(b)
        matrix __Vybar_`s' = e(V)

        // Under Lee, a negative rho implies POSITIVE selection bias, because
        // the correction term is strictly positive (Lee 2007). Confirm the
        // sign convention of your installed version of selmlog before
        // interpreting this quantity.
        di as txt "  theta = " as res %9.6f __th_`s'[1,1] ///
           as txt "   hbar = " as res %9.6f __hbar_`s'[1,1] ///
           as txt "   rho = "  as res %9.6f `rho_`s''
    }

    // ======================================================================
    // 3. DECOMPOSITION
    // ======================================================================

    mata: __ob_lee_decomp()

    tempname B S
    matrix `B' = __res_b
    matrix `S' = __res_se
    local cnames "Difference Endowments Returns Selectivity WageOffers"
    matrix colnames `B' = `cnames'
    matrix colnames `S' = `cnames'

    // NOTE: -ereturn post- and -ereturn matrix- ERASE the source matrix.
    // Anything needed afterwards must be read back from e().
    ereturn post `B', esample(`touse') depname("`depvar'")
    ereturn matrix se_delta   = `S'
    ereturn scalar closure    = __res_chk    // diff - (endow+retn+selec)
    ereturn scalar closure2   = __res_chk2   // offers - (diff - selec)
    ereturn scalar N_ref      = `nobs_m'
    ereturn scalar N_other    = `nobs_f'
    ereturn scalar theta_ref  = __th_m[1,1]
    ereturn scalar theta_oth  = __th_f[1,1]
    ereturn scalar hbar_ref   = __hbar_m[1,1]
    ereturn scalar hbar_oth   = __hbar_f[1,1]
    ereturn scalar rho_ref    = `rho_m'
    ereturn scalar rho_oth    = `rho_f'
    ereturn scalar var_diff   = __res_se[1,1]^2   // input to the Welch test
    ereturn local  ref        "`gm'"
    ereturn local  other      "`gf'"
    ereturn local  secondstep = cond("`nowls'" == "", "wls", "ols")
    ereturn local  vcetype    = cond("`vrobust'" == "", "conventional", "robust")
    ereturn local  cmd        "oaxaca_twofold_lee_v2"

    if "`nodetail'" == "" {
        tempname D
        matrix `D' = __res_det
        matrix rownames `D' = `x_ols' _cons
        matrix colnames `D' = Endowments Returns
        ereturn matrix detailed = `D'
    }

    // ======================================================================
    // 4. DISPLAY (read back from e(); the tempnames above no longer exist)
    // ======================================================================

    tempname BD SD
    matrix `BD' = e(b)
    matrix `SD' = e(se_delta)

    di as txt _n "{hline 74}"
    di as txt "Decomposition: Neuman-Oaxaca (2004) noncommittal, Lee correction"
    di as txt "{hline 74}"
    di as txt %14s "Component" %14s "Estimate" %14s "Std. err." ///
              %15s "CI lower" %15s "CI upper"
    di as txt "{hline 74}"
    forvalues j = 1/5 {
        local nm : word `j' of `cnames'
        local b_j = `BD'[1, `j']
        local s_j = `SD'[1, `j']
        if `j' == 5 di as txt "{hline 74}"
        di as txt %14s "`nm'" as res %14.6f `b_j' %14.6f `s_j' ///
           %15.6f (`b_j' - 1.96*`s_j') %15.6f (`b_j' + 1.96*`s_j')
    }
    di as txt "{hline 74}"
    di as txt "Components 1-4 add up to the observed gap."
    di as txt "WageOffers = Xbar_ref'b_ref - Xbar_other'b_other, the gap that"
    di as txt "would prevail in the absence of differential selection."
    di as txt "Closure 1 (accounting identity) : " as res %18.16f e(closure)
    di as txt "Closure 2 (WageOffers routes)   : " as res %18.16f e(closure2)

    //----------------------------------------------------------------------
    // 4.1 Stability diagnostic
    //----------------------------------------------------------------------
    local flagged 0
    forvalues j = 2/4 {
        local nm : word `j' of `cnames'
        local b_j = abs(`BD'[1, `j'])
        local s_j = `SD'[1, `j']
        if (`b_j' > 0) & (`s_j' > 5 * `b_j') {
            if `flagged' == 0 {
                di as txt ""
                di as err "UNSTABLE: the delta-method standard error exceeds"
                di as err "five times the estimate in absolute value for:"
                local flagged 1
            }
            di as err "    `nm' (se/|b| = " %6.1f (`s_j'/`b_j') ")"
        }
    }
    if `flagged' {
        di as err ""
        di as err "This is the signature of an index whose influence function"
        di as err "has heavy tails: GE(-1) involves y^(-1), its variance"
        di as err "depends on E[y^(-2)] and explodes when wages are low."
        di as err "Do not report these intervals. Use -bootstrap-; if the"
        di as err "width survives, it is a result about precision."
    }

    //----------------------------------------------------------------------
    // 4.2 Accounting-closure diagnostic
    //----------------------------------------------------------------------
    if abs(e(closure)) > 1e-8 {
        di as txt ""
        if "`nowls'" == "" {
            di as err "The closure is not exact: the second stage is weighted"
            di as err "and the means are not. Re-run with -nowls-."
        }
        else {
            di as err "The closure is not exact under OLS. Either a block of"
            di as err "coefficients is missing, or the sample does not match"
            di as err "the second stage. Compare N_ref/N_other with the"
            di as err "regression."
        }
    }
    if abs(e(closure2)) > 1e-8 {
        di as err ""
        di as err "The two routes to WageOffers disagree. This points to a"
        di as err "problem in the extraction of theta or of hbar."
    }

    if "`vrobust'" == "" {
        di as txt _n "Note: standard errors use the nonrobust e(V) of selmlog."
        di as txt    "If the uncorrected branch of your analysis uses robust"
        di as txt    "errors, the two are not comparable. Add -vrobust-."
    }
    di as txt _n "All standard errors treat the first-stage probabilities as"
    di as txt    "known and are a LOWER BOUND. Publishable inference should"
    di as txt    "come from -bootstrap, reps(500)-."

    capture matrix drop __b_m __b_f __V_bm __V_bf __th_m __th_f
    capture matrix drop __V_thm __V_thf __xbar_m __xbar_f
    capture matrix drop __Vxbar_m __Vxbar_f __hbar_m __hbar_f
    capture matrix drop __Vhbar_m __Vhbar_f __ybar_m __ybar_f
    capture matrix drop __Vybar_m __Vybar_f
    capture matrix drop __res_b __res_se __res_det
    capture scalar drop __res_chk __res_chk2
end


// =========================================================================
// MATA: a single function; inputs and outputs through named Stata matrices
// =========================================================================

mata:

void __ob_lee_decomp()
{
    real colvector bm, bf, xm, xf, thm, thf, hm, hf
    real colvector dx, db, endow_d, retn_d
    real matrix Vbm, Vbf, Vthm, Vthf, Vxm, Vxf, Vhm, Vhf
    real scalar endow, retn, selec, diff, offers, chk, chk2
    real scalar vendow, vretn, vselec, vdiff, voffers

    // --- input: suffix m = reference, f = comparison ---------------------
    bm  = st_matrix("__b_m")'          // beta of the reference group
    bf  = st_matrix("__b_f")'
    Vbm = st_matrix("__V_bm")
    Vbf = st_matrix("__V_bf")

    thm  = st_matrix("__th_m")'        // theta (coefficient on the Lee term)
    thf  = st_matrix("__th_f")'
    Vthm = st_matrix("__V_thm")
    Vthf = st_matrix("__V_thf")

    xm  = st_matrix("__xbar_m")'       // covariate means
    xf  = st_matrix("__xbar_f")'
    Vxm = st_matrix("__Vxbar_m")
    Vxf = st_matrix("__Vxbar_f")

    hm  = st_matrix("__hbar_m")'       // means of the correction term
    hf  = st_matrix("__hbar_f")'
    Vhm = st_matrix("__Vhbar_m")
    Vhf = st_matrix("__Vhbar_f")

    dx = xm - xf                       // endowment differential
    db = bm - bf                       // returns differential

    // --- components: Neuman and Oaxaca (2004), noncommittal scheme -------
    endow = dx' * bm                   // endowments, valued at beta_ref
    retn  = xf' * db                   // returns, evaluated at Xbar_other
    selec = thm' * hm - thf' * hf      // selectivity
    diff  = st_matrix("__ybar_m")[1,1] - st_matrix("__ybar_f")[1,1]

    // Accounting closure: zero to machine precision under OLS
    chk   = diff - (endow + retn + selec)

    // --- wage-offer gap --------------------------------------------------
    // By the accounting identity,
    //   offers = diff - selec = endow + retn = Xbar_r'b_r - Xbar_o'b_o.
    // It is computed from the last expression because that way it depends
    // ONLY on (beta, Xbar) of both groups, and its variance has an exact form
    // within the delta framework, with no discarded covariances.
    offers = xm' * bm - xf' * bf

    // Second correctness check: the two routes must agree
    chk2 = offers - (diff - selec)

    // --- detailed components by covariate --------------------------------
    endow_d = dx :* bm
    retn_d  = xf :* db

    // --- delta-method variances ------------------------------------------
    // Difference: the two means come from independent samples.
    vdiff = st_matrix("__Vybar_m")[1,1] + st_matrix("__Vybar_f")[1,1]

    // Endowments and returns: the first-order approximation of Jann (2008),
    // which treats beta and Xbar as independent of each other and across
    // groups.
    vendow = dx' * Vbm * dx + bm' * Vxm * bm + bm' * Vxf * bm
    vretn  = xf' * (Vbm + Vbf) * xf + db' * Vxf * db

    // Selectivity: the same logic applied to theta and hbar.
    vselec = hm' * Vthm * hm + thm' * Vhm * thm +
             hf' * Vthf * hf + thf' * Vhf * thf

    // Offers: exact within the delta framework.
    voffers = xm' * Vbm * xm + bm' * Vxm * bm +
              xf' * Vbf * xf + bf' * Vxf * bf

    // --- output ------------------------------------------------------------
    st_matrix("__res_b",  (diff, endow, retn, selec, offers))
    st_matrix("__res_se", (sqrt(vdiff), sqrt(vendow), sqrt(vretn),
                           sqrt(vselec), sqrt(voffers)))
    st_matrix("__res_det", (endow_d, retn_d))
    st_numscalar("__res_chk",  chk)
    st_numscalar("__res_chk2", chk2)
}

end

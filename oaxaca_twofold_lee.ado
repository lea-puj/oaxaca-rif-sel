*! oaxaca_twofold_lee_v2 v2.9  -- 2026-08-23
*! ============================================================================
*! DESCOMPOSICION DE OAXACA-BLINDER EN DOS COMPONENTES CON CORRECCION DE
*! SESGO DE SELECCION POLICOTOMICA DE LEE (1983), REPORTADA BAJO EL ESQUEMA
*! NO COMPROMETIDO DE NEUMAN Y OAXACA (2004, ecuacion 11).
*!
*! Autores del proyecto: Mauricio Quinones, Sergio Barona
*! Documentacion formal completa: nota_tecnica_ado.tex
*!
*! ----------------------------------------------------------------------------
*! QUE REPORTA
*! ----------------------------------------------------------------------------
*!   Difference   brecha observada, Ybar(ref) - Ybar(otro)
*!   Endowments   (Xbar_ref - Xbar_otro)' beta_ref        [dotaciones]
*!   Returns      Xbar_otro' (beta_ref - beta_otro)       [retornos]
*!   Selectivity  theta_ref*hbar_ref - theta_otro*hbar_otro
*!   WageOffers   Xbar_ref'beta_ref - Xbar_otro'beta_otro
*!
*!   Los primeros cuatro son la ecuacion (11) de Neuman-Oaxaca y suman
*!   exactamente la brecha observada. WageOffers es el estimando de su
*!   ecuacion (7): la brecha que prevaleceria en ausencia de seleccion
*!   diferencial (Ben Yahmed 2018, World Development, Tabla 9, panel 3.1).
*!
*! ----------------------------------------------------------------------------
*! HISTORIAL DE VERSIONES
*! ----------------------------------------------------------------------------
*! v2.9  CORRECCION: Stata antepone "o." al nombre de los coeficientes omitidos
*!       por colinealidad, de modo que e(b) puede contener "o.sector_5" en vez
*!       de "sector_5". La busqueda por nombre exacto fallaba con rc=111. Los
*!       sectores colineales varian entre celdas (en informales se omiten
*!       sector_5 y sector_16, en formales no), asi que el fallo se producia
*!       solo en algunas celdas. Ahora se limpia cualquier prefijo de factor
*!       ("o.", "b.", "1b.", "2o.") antes de comparar. Un coeficiente omitido
*!       vale cero y aporta cero a ambos componentes, que es lo correcto; su
*!       inclusion solo mantiene alineadas las dimensiones de beta y Xbar.
*!
*!       NOTA PARA EL APENDICE: como el conjunto de sectores omitidos difiere
*!       entre celdas, los componentes DETALLADOS por covariable no son
*!       directamente comparables entre sectores. Los agregados si lo son,
*!       porque el termino explicado suma sobre todas las dummies y es
*!       invariante a cual se omita.
*! v2.8  Varianza exacta de WageOffers y segundo cierre e(closure2).
*! v2.7  Corrige el error r(111) por leer tempnames tras -ereturn post-.
*! v2.6  Error estandar real para Difference; vrobust; diagnosticos.
*! v2.5  Muestra reconstruida con markout.
*! v2.4  Opcion nowls: cierre contable exacto.
*! v2.3  Terminos de correccion localizados por el prefijo _m.
*! v2.2  La primera etapa conserva los tres estados.
*!
*! ----------------------------------------------------------------------------
*! SINTAXIS
*! ----------------------------------------------------------------------------
*!   oaxaca_twofold_lee_v2 depvar indepvars [if] [in] [aw pw fw], ///
*!       by(grupo) select(zvarlist) [ref(#) status(nombre) genpfx(nombre)
*!        innerboot(#) nowls vrobust nodetail]
*!
*!   depvar   DEBE estar en missing fuera del sector de interes. NO pasar
*!            -if status_var==#-: la primera etapa necesita los tres estados.
*!   nowls    segunda etapa por MCO. RECOMENDADO: hace exacto el cierre.
*!   vrobust  varianza robusta a heterocedasticidad; implica nowls.
*! ============================================================================

program define oaxaca_twofold_lee_v2, eclass
    version 16

    // ======================================================================
    // 0. SINTAXIS Y OPCIONES
    // ======================================================================

    syntax varlist(min=2) [if] [in] [aweight pweight fweight],           ///
        BY(varname) SELect(varlist)                                      ///
        [ REF(string) STATus(name) GENpfx(name)                          ///
          INNERboot(integer 0) NOWLS VROBust NODEtail ]

    gettoken depvar indepvars : varlist
    local x_ols "`indepvars'"

    if "`status'" == "" local status status_var
    if "`genpfx'" == "" local genpfx _h

    // vrobust toma b y V de un reajuste por MCO, de modo que la segunda etapa
    // debe ser MCO para que los coeficientes reportados sean los de selmlog.
    if "`vrobust'" != "" & "`nowls'" == "" {
        di as txt "Nota: vrobust implica segunda etapa por MCO; se asume nowls."
        local nowls "nowls"
    }
    local wlsopt = cond("`nowls'" == "", "wls", "")

    // touse NO debe excluir observaciones por tener `depvar' en missing: esas
    // observaciones identifican la primera etapa. El recorte sectorial viaja
    // por la dependiente, no por la muestra.
    marksample touse, novarlist
    markout `touse' `by' `select' `status'

    // ======================================================================
    // 1. GRUPOS Y ESTRUCTURA DE REFERENCIA
    // ======================================================================

    qui levelsof `by' if `touse', local(glevels)
    local ng : word count `glevels'
    if `ng' != 2 {
        di as err "`by' debe tomar exactamente dos valores en la muestra."
        exit 198
    }
    local lev1 : word 1 of `glevels'
    local lev2 : word 2 of `glevels'

    if "`ref'" == "" {
        local gm `lev2'
        di as txt "Nota: no se especifico ref(); se usa `by' == `gm'."
    }
    else {
        local gm `ref'
        if (`gm' != `lev1') & (`gm' != `lev2') {
            di as err "ref(`gm') no es un valor de `by'."
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
    di as txt _n "Estructura de referencia : `by' == `gm' [`gm_lbl']"
    di as txt    "Grupo de comparacion     : `by' == `gf' [`gf_lbl']"
    di as txt    "Brecha reportada         : Ybar(`gm_lbl') - Ybar(`gf_lbl')"
    di as txt    "Segunda etapa            : " ///
                 cond("`nowls'" == "", "MCP (default de selmlog)", "MCO")
    di as txt    "Varianza                 : " ///
                 cond("`vrobust'" == "", "e(V) de selmlog, no robusta", ///
                                         "robusta a heterocedasticidad") _n

    local bootopt ""
    if `innerboot' > 0 local bootopt "bootstrap(`innerboot')"

    tempvar one
    qui gen byte `one' = 1

    // ======================================================================
    // 2. ESTIMACION POR GRUPO
    // ======================================================================
    // Sufijo "m" = grupo de referencia, "f" = grupo de comparacion.
    // ======================================================================

    foreach s in m f {

        if "`s'" == "m"  local gval `gm'
        else             local gval `gf'

        //------------------------------------------------------------------
        // 2.1 Verificacion: la primera etapa necesita los tres estados
        //------------------------------------------------------------------
        qui levelsof `status' if `touse' & `by' == `gval', local(_stlev)
        local nst : word count `_stlev'
        if `nst' < 3 {
            di as err _n "`status' toma solo `nst' valor(es) en esta muestra."
            di as err "La correccion de Lee necesita los tres estados en la"
            di as err "PRIMERA ETAPA. No pase -if `status'==#-. Restrinja el"
            di as err "sector a traves de la variable dependiente:"
            di as err "    gen dep_formal = dep_var if `status'==1"
            exit 148
        }

        qui count if `touse' & `by' == `gval' & !missing(`depvar')
        if r(N) < 30 {
            di as err "Solo `r(N)' observaciones con `depvar' no faltante" ///
                      " para `by' == `gval'."
            exit 2001
        }

        //------------------------------------------------------------------
        // 2.2 Estimacion en dos etapas (selmlog)
        //------------------------------------------------------------------
        // La advertencia de selmlog sobre observaciones con valores faltantes
        // en la ecuacion principal es el comportamiento buscado: es lo que
        // hace identificable la correccion.
        //------------------------------------------------------------------
        capture drop `genpfx'*
        di as txt "Estimando para `by' == `gval' ..."

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
        // 2.3 Localizacion de los terminos de correccion
        //------------------------------------------------------------------
        // selmlog nombra el COEFICIENTE como _m<j> dentro de la ecuacion
        // salarial, y coloca Sigma<j> y rho<j> en una ecuacion aparte llamada
        // "Anciliary". La VARIABLE correspondiente en la base se llama
        // `genpfx'<j>, por la opcion gen(). Se enlazan por el sufijo numerico:
        // el coeficiente da theta_j, la variable da hbar_j.
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
                    di as err "Se encontro el coeficiente `nm' pero la" ///
                              " variable `genpfx'`jsuf' no existe."
                    di as err "Revise el prefijo pasado en genpfx()."
                    exit 111
                }
                local hvars "`hvars' `genpfx'`jsuf'"
            }
        }
        if "`hvars'" == "" {
            di as err _n "No se hallaron terminos de correccion (_m*) en e(b)."
            di as err "Nombres de columna: `cn'"
            exit 111
        }
        di as txt "  termino(s) de correccion: `hvars'"

        //------------------------------------------------------------------
        // 2.4 Muestra de la ecuacion salarial
        //------------------------------------------------------------------
        // No se puede usar e(sample) de selmlog: bajo -nowls- viene vacia.
        // Se reconstruye reproduciendo la eliminacion por lista del segundo
        // paso. El conteo impreso debe coincidir con el "Number of obs" de la
        // regresion de segunda etapa.
        //------------------------------------------------------------------
        tempvar es`s'
        qui gen byte `es`s'' = (`touse' & `by' == `gval' & !missing(`depvar'))
        markout `es`s'' `x_ols' `hvars'
        qui count if `es`s''
        local nobs_`s' = r(N)
        di as txt "  observaciones de la ecuacion salarial: `nobs_`s''"
        if `nobs_`s'' == 0 {
            di as err "Muestra vacia para `by' == `gval'."
            exit 2000
        }

        //------------------------------------------------------------------
        // 2.5 Reajuste robusto (opcional)
        //------------------------------------------------------------------
        // Bajo nowls la segunda etapa es MCO de depvar sobre x_ols y los
        // terminos de correccion. Reajustar con -regress, robust- reproduce
        // exactamente los coeficientes de selmlog y solo reemplaza la matriz
        // de varianzas. Importa porque la rama sin correccion (oaxaca_rif)
        // usa errores robustos: sin esto, los dos paneles de cada figura no
        // son comparables en inferencia.
        //------------------------------------------------------------------
        if "`vrobust'" != "" {
            qui regress `depvar' `x_ols' `hvars' if `es`s'', robust
            matrix `bfull' = e(b)
            matrix `Vfull' = e(V)
            local cn : colnames `bfull'
            local ce : coleq    `bfull'
            local k = colsof(`bfull')
            local hnames "`hvars'"      // regress los nombra por variable
        }
        else {
            local hnames ""             // selmlog los nombra _m<j>
        }

        //------------------------------------------------------------------
        // 2.6 Posiciones de las covariables
        //------------------------------------------------------------------
        // ATENCION AL PREFIJO "o.": Stata lo antepone al nombre de los
        // coeficientes omitidos por colinealidad, de modo que e(b) puede
        // contener "o.sector_5" en lugar de "sector_5". Comparar por nombre
        // exacto fallaba con rc=111, y como los sectores colineales varian
        // entre celdas (en informales se omiten sector_5 y sector_16), el
        // fallo aparecia solo en algunas celdas.
        //
        // regexr limpia cualquier prefijo de factor: "o.", "b.", "1b.", "2o.".
        // Un coeficiente omitido vale cero y aporta cero tanto a dotaciones
        // como a retornos, que es lo correcto; incluirlo solo mantiene
        // alineadas las dimensiones de beta y Xbar. El cierre contable lo
        // confirma.
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
                di as err "No se hallo el coeficiente `v'."
                di as err "Columnas: `cn'"
                exit 111
            }
        }
        if `nomit' > 0 {
            di as txt "  covariables omitidas por colinealidad: `nomit'" ///
                      " (coeficiente cero, aportan cero a la descomposicion)"
        }

        //------------------------------------------------------------------
        // 2.7 Posiciones de los terminos de correccion
        //------------------------------------------------------------------
        local tpos ""
        if "`hnames'" != "" {
            // Ruta con vrobust: buscar por nombre de variable
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
                    di as err "No se hallo el termino `v' tras el reajuste."
                    di as err "Columnas: `cn'"
                    exit 111
                }
            }
        }
        else {
            // Ruta sin vrobust: buscar por el prefijo _m del coeficiente
            forvalues c = 1/`k' {
                local nm : word `c' of `cn'
                local eq : word `c' of `ce'
                if inlist(lower("`eq'"), "anciliary", "ancillary") continue
                local nm_lim = regexr("`nm'", "^[0-9]*[bo]?\.", "")
                if substr("`nm_lim'", 1, 2) == "_m" local tpos "`tpos'`c' "
            }
        }

        //------------------------------------------------------------------
        // 2.8 Extraccion de submatrices
        //------------------------------------------------------------------
        // Se extraen dos bloques de e(b) y e(V): covariables (beta, V_beta) y
        // terminos de correccion (theta, V_theta). Se construyen elemento por
        // elemento porque Stata no indexa matrices por vector de posiciones.
        //
        // Se descarta la covarianza entre beta y theta (bloques fuera de la
        // diagonal). Es inocuo: ningun componente individual depende de ambos.
        //------------------------------------------------------------------
        tempname bsub Vsub tsub Vtsub

        // --- beta y su varianza ---
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

        // --- theta y su varianza ---
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
        // 2.9 Medias sobre la muestra de la ecuacion salarial
        //------------------------------------------------------------------
        // -mean- devuelve en e(b) el vector de medias y en e(V) su matriz de
        // varianzas-covarianzas, que es exactamente V(Xbar), respetando los
        // pesos. La media de la dependiente aporta Ybar y su varianza, que es
        // lo que necesita el componente Difference.
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

        // Bajo Lee, un rho negativo implica sesgo de seleccion POSITIVO,
        // porque el termino de correccion es estrictamente positivo (Lee
        // 2007). Verifique la convencion de signo de selmlog antes de
        // interpretarlo en el texto.
        di as txt "  theta = " as res %9.6f __th_`s'[1,1] ///
           as txt "   hbar = " as res %9.6f __hbar_`s'[1,1] ///
           as txt "   rho = "  as res %9.6f `rho_`s''
    }

    // ======================================================================
    // 3. DESCOMPOSICION
    // ======================================================================

    mata: __ob_lee_decomp()

    tempname B S
    matrix `B' = __res_b
    matrix `S' = __res_se
    local cnames "Difference Endowments Returns Selectivity WageOffers"
    matrix colnames `B' = `cnames'
    matrix colnames `S' = `cnames'

    // ATENCION: -ereturn post- y -ereturn matrix- BORRAN la matriz de origen.
    // Todo lo que se necesite despues hay que releerlo desde e().
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
    ereturn scalar var_diff   = __res_se[1,1]^2   // insumo del test de Welch
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
    // 4. DESPLIEGUE (releyendo desde e(); los tempnames de arriba ya no existen)
    // ======================================================================

    tempname BD SD
    matrix `BD' = e(b)
    matrix `SD' = e(se_delta)

    di as txt _n "{hline 74}"
    di as txt "Descomposicion: Neuman-Oaxaca (2004) ec. (11), correccion de Lee"
    di as txt "{hline 74}"
    di as txt %14s "Componente" %14s "Estimador" %14s "Error est." ///
              %15s "IC inferior" %15s "IC superior"
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
    di as txt "Los componentes 1-4 son la ec. (11) y suman la brecha observada."
    di as txt "WageOffers = Xbar_ref'b_ref - Xbar_otro'b_otro, el estimando de"
    di as txt "la ec. (7) (Ben Yahmed 2018, Tabla 9, panel 3.1)."
    di as txt "Cierre 1 (identidad contable) : " as res %18.16f e(closure)
    di as txt "Cierre 2 (rutas de WageOffers): " as res %18.16f e(closure2)

    //----------------------------------------------------------------------
    // 4.1 Diagnostico de estabilidad
    //----------------------------------------------------------------------
    local flagged 0
    forvalues j = 2/4 {
        local nm : word `j' of `cnames'
        local b_j = abs(`BD'[1, `j'])
        local s_j = `SD'[1, `j']
        if (`b_j' > 0) & (`s_j' > 5 * `b_j') {
            if `flagged' == 0 {
                di as txt ""
                di as err "INESTABLE: el error estandar delta-method supera"
                di as err "cinco veces el estimador en valor absoluto para:"
                local flagged 1
            }
            di as err "    `nm' (se/|b| = " %6.1f (`s_j'/`b_j') ")"
        }
    }
    if `flagged' {
        di as err ""
        di as err "Es la firma de un indice cuya funcion de influencia tiene"
        di as err "colas pesadas: GE(-1) involucra y^(-1), su varianza depende"
        di as err "de E[y^(-2)] y explota cuando los salarios son bajos."
        di as err "No reporte estos intervalos. Use -bootstrap-; si la"
        di as err "amplitud sobrevive, es un resultado sobre la precision."
    }

    //----------------------------------------------------------------------
    // 4.2 Diagnostico de cierre contable
    //----------------------------------------------------------------------
    if abs(e(closure)) > 1e-8 {
        di as txt ""
        if "`nowls'" == "" {
            di as err "El cierre no es exacto: la segunda etapa es ponderada y"
            di as err "las medias no lo son. Vuelva a correr con -nowls-."
        }
        else {
            di as err "El cierre no es exacto bajo MCO. O falta un bloque de"
            di as err "coeficientes, o la muestra no coincide con la segunda"
            di as err "etapa. Compare N_ref/N_other con la regresion."
        }
    }
    if abs(e(closure2)) > 1e-8 {
        di as err ""
        di as err "Las dos rutas de calculo de WageOffers no coinciden."
        di as err "Indica un problema en la extraccion de theta o de hbar."
    }

    if "`vrobust'" == "" {
        di as txt _n "Nota: los errores usan e(V) no robusta de selmlog. La"
        di as txt    "rama sin correccion (oaxaca_rif) usa errores robustos,"
        di as txt    "de modo que ambas no son comparables. Agregue -vrobust-."
    }
    di as txt _n "Todos los errores estandar tratan las probabilidades de la"
    di as txt    "primera etapa como conocidas y son una COTA INFERIOR. La"
    di as txt    "inferencia publicable debe venir de -bootstrap, reps(500)-."

    capture matrix drop __b_m __b_f __V_bm __V_bf __th_m __th_f
    capture matrix drop __V_thm __V_thf __xbar_m __xbar_f
    capture matrix drop __Vxbar_m __Vxbar_f __hbar_m __hbar_f
    capture matrix drop __Vhbar_m __Vhbar_f __ybar_m __ybar_f
    capture matrix drop __Vybar_m __Vybar_f
    capture matrix drop __res_b __res_se __res_det
    capture scalar drop __res_chk __res_chk2
end


// =========================================================================
// MATA: una sola funcion; entradas y salidas por matrices de Stata con nombre
// =========================================================================

mata:

void __ob_lee_decomp()
{
    real colvector bm, bf, xm, xf, thm, thf, hm, hf
    real colvector dx, db, endow_d, retn_d
    real matrix Vbm, Vbf, Vthm, Vthf, Vxm, Vxf, Vhm, Vhf
    real scalar endow, retn, selec, diff, offers, chk, chk2
    real scalar vendow, vretn, vselec, vdiff, voffers

    // --- lectura: sufijo m = referencia, f = comparacion -----------------
    bm  = st_matrix("__b_m")'          // beta del grupo de referencia
    bf  = st_matrix("__b_f")'
    Vbm = st_matrix("__V_bm")
    Vbf = st_matrix("__V_bf")

    thm  = st_matrix("__th_m")'        // theta (coef. del termino de Lee)
    thf  = st_matrix("__th_f")'
    Vthm = st_matrix("__V_thm")
    Vthf = st_matrix("__V_thf")

    xm  = st_matrix("__xbar_m")'       // medias de covariables
    xf  = st_matrix("__xbar_f")'
    Vxm = st_matrix("__Vxbar_m")
    Vxf = st_matrix("__Vxbar_f")

    hm  = st_matrix("__hbar_m")'       // medias del termino de correccion
    hf  = st_matrix("__hbar_f")'
    Vhm = st_matrix("__Vhbar_m")
    Vhf = st_matrix("__Vhbar_f")

    dx = xm - xf                       // diferencial de dotaciones
    db = bm - bf                       // diferencial de retornos

    // --- componentes: Neuman-Oaxaca (2004), ecuacion (11) ----------------
    endow = dx' * bm                   // dotaciones, valoradas a beta_ref
    retn  = xf' * db                   // retornos, evaluados en Xbar_otro
    selec = thm' * hm - thf' * hf      // selectividad
    diff  = st_matrix("__ybar_m")[1,1] - st_matrix("__ybar_f")[1,1]

    // Cierre contable: cero a precision de maquina bajo MCO
    chk   = diff - (endow + retn + selec)

    // --- brecha en ofertas salariales: ecuacion (7) -----------------------
    // Por la identidad contable,
    //   offers = diff - selec = endow + retn = Xbar_r'b_r - Xbar_o'b_o.
    // Se calcula por la ultima expresion porque asi depende SOLO de
    // (beta, Xbar) de ambos grupos, y su varianza tiene forma exacta dentro
    // del marco delta, sin covarianzas descartadas.
    offers = xm' * bm - xf' * bf

    // Segunda prueba de correctitud: las dos rutas deben coincidir
    chk2 = offers - (diff - selec)

    // --- componentes detallados por covariable ---------------------------
    endow_d = dx :* bm
    retn_d  = xf :* db

    // --- varianzas por el metodo delta -----------------------------------
    // Difference: las dos medias provienen de muestras independientes.
    vdiff = st_matrix("__Vybar_m")[1,1] + st_matrix("__Vybar_f")[1,1]

    // Dotaciones y retornos: aproximacion de primer orden de Jann (2008),
    // que trata beta y Xbar como independientes entre si y entre grupos.
    vendow = dx' * Vbm * dx + bm' * Vxm * bm + bm' * Vxf * bm
    vretn  = xf' * (Vbm + Vbf) * xf + db' * Vxf * db

    // Selectividad: misma logica sobre theta y hbar.
    vselec = hm' * Vthm * hm + thm' * Vhm * thm +
             hf' * Vthf * hf + thf' * Vhf * thf

    // Ofertas: exacta dentro del marco delta.
    voffers = xm' * Vbm * xm + bm' * Vxm * bm +
              xf' * Vbf * xf + bf' * Vxf * bf

    // --- salida -----------------------------------------------------------
    st_matrix("__res_b",  (diff, endow, retn, selec, offers))
    st_matrix("__res_se", (sqrt(vdiff), sqrt(vendow), sqrt(vretn),
                           sqrt(vselec), sqrt(voffers)))
    st_matrix("__res_det", (endow_d, retn_d))
    st_numscalar("__res_chk",  chk)
    st_numscalar("__res_chk2", chk2)
}

end

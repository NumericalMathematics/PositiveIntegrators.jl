var documenterSearchIndex = {"docs":
[{"location":"api_reference/#PositiveIntegrators.jl-API","page":"API reference","title":"PositiveIntegrators.jl API","text":"","category":"section"},{"location":"api_reference/","page":"API reference","title":"API reference","text":"CurrentModule = PositiveIntegrators","category":"page"},{"location":"api_reference/","page":"API reference","title":"API reference","text":"Modules = [PositiveIntegrators]","category":"page"},{"location":"api_reference/#PositiveIntegrators.prob_pds_bertolazzi","page":"API reference","title":"PositiveIntegrators.prob_pds_bertolazzi","text":"prob_pds_bertolazzi\n\nPositive and conservative autonomous nonlinear PDS\n\nbeginaligned\nmathbfu=beginpmatrix2 -1 -1-1 2 -1-1 -1 2endpmatrixbeginpmatrix5u_2u_3(10^-2 + (u_2u_3)^2) + u_2u_3(10^-16 + u_2u_3(10^-8 + u_2u_3))\n10u_1u_3^2\n01(u_3 - u_2 - 25)^2u_1u_2endpmatrix\nendaligned\n\nwith initial value mathbfu_0 = (00 10 20)^T and time domain (00 10). There is one independent linear invariant, e.g. u_1+u_2+u_3 = 30.\n\nReferences\n\nEnrico Bertolazzi. \"Positive and conservative schemes for mass action kinetics.\" Computers and Mathematics with Applications 32 (1996): 29-43. DOI: 10.1016/0898-1221(96)00142-3\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.prob_pds_brusselator","page":"API reference","title":"PositiveIntegrators.prob_pds_brusselator","text":"prob_pds_brusselator\n\nPositive and conservative autonomous nonlinear PDS\n\nbeginaligned\nu_1 = -u_1\nu_2 = -u_2u_5\nu_3 = u_2u_5\nu_4 = u_5\nu_5 = u_1 - u_2u_5 + u_5^2u_6 - u_5\nu_6 = u_2u_5 - u_5^2u_6\nendaligned\n\nwith initial value mathbfu_0 = (100 100 00 00 01 01)^T and time domain (00 200). There are two independent linear invariants, e.g. u_1+u_4+u_5+u_6 = 102 and u_2+u_3 = 100.\n\nReferences\n\nLuca Bonaventura,  and Alessandro Della Rocca. \"Unconditionally Strong Stability Preserving Extensions of the TR-BDF2 Method.\" Journal of Scientific Computing 70 (2017): 859 - 895. DOI: 10.1007/s10915-016-0267-9\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.prob_pds_linmod","page":"API reference","title":"PositiveIntegrators.prob_pds_linmod","text":"prob_pds_linmod\n\nPositive and conservative autonomous linear PDS\n\nbeginaligned\nu_1 = u_2 - 5u_1\nu_2 = 5u_1 - u_2\nendaligned\n\nwith initial value mathbfu_0 = (09 01)^T and time domain (00 20). There is one independent linear invariant, e.g. u_1+u_2 = 1.\n\nReferences\n\nHans Burchard, Eric Deleersnijder, and Andreas Meister. \"A high-order conservative Patankar-type discretisation for stiff systems of production-destruction equations.\" Applied Numerical Mathematics 47.1 (2003): 1-30. DOI: 10.1016/S0168-9274(03)00101-6\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.prob_pds_linmod_inplace","page":"API reference","title":"PositiveIntegrators.prob_pds_linmod_inplace","text":"prob_pds_linmod_inplace\n\nSame as prob_pds_linmod but with in-place computation.\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.prob_pds_nonlinmod","page":"API reference","title":"PositiveIntegrators.prob_pds_nonlinmod","text":"prob_pds_nonlinmod\n\nPositive and conservative autonomous nonlinear PDS\n\nbeginaligned\nu_1 = -fracu_1u_2u_1 + 10\nu_2 = fracu_1u_2u_1 + 10 - 03u_2\nu_3 = 03 u_2\nendaligned\n\nwith initial value mathbfu_0 = (998 001 001)^T and time domain (00 300). There is one independent linear invariant, e.g. u_1+u_2+u_3 = 100.\n\nReferences\n\nHans Burchard, Eric Deleersnijder, and Andreas Meister. \"A high-order conservative Patankar-type discretisation for stiff systems of production-destruction equations.\" Applied Numerical Mathematics 47.1 (2003): 1-30. DOI: 10.1016/S0168-9274(03)00101-6\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.prob_pds_npzd","page":"API reference","title":"PositiveIntegrators.prob_pds_npzd","text":"prob_pds_npzd\n\nPositive and conservative autonomous nonlinear PDS\n\nbeginaligned\nu_1 = 001u_2 + 001u_3 + 0003u_4 - fracu_1u_2001 + u_1\nu_2 = fracu_1u_2001 + u_1- 001u_2 - 05( 1 - e^-121u_2^2)u_3 - 005u_2\nu_3 = 05(1 - e^-121u_2^2)u_3 - 001u_3 - 002u_3\nu_4 = 005u_2 + 002u_3 - 0003u_4\nendaligned\n\nwith initial value mathbfu_0 = (80 20 10 40)^T and time domain (00 100). There is one independent linear invariant, e.g. u_1+u_2+u_3+u_4 = 150.\n\nReferences\n\nHans Burchard, Eric Deleersnijder, and Andreas Meister. \"Application of modified Patankar schemes to stiff biogeochemical models for the water column.\" Ocean Dynamics 55 (2005): 326-337. DOI: 10.1007/s10236-005-0001-x\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.prob_pds_robertson","page":"API reference","title":"PositiveIntegrators.prob_pds_robertson","text":"prob_pds_robertson\n\nPositive and conservative autonomous nonlinear PDS\n\nbeginaligned\nu_1 = -004u_1+10^4 u_2u_3\nu_2 =  004u_1-10^4 u_2u_3-310^7 u_2^2\nu_3 = 310^7 u_2^2\nendaligned\n\nwith initial value mathbfu_0 = (10 00 00)^T and time domain (00 10^11). There is one independent linear invariant, e.g. u_1+u_2+u_3 = 10.\n\nReferences\n\nErnst Hairer, Gerd Wanner. \"Solving Ordinary Differential Equations II - Stiff and Differential-Algebraic Problems.\" 2nd Edition, Springer (2002): Section IV.1.\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.prob_pds_sir","page":"API reference","title":"PositiveIntegrators.prob_pds_sir","text":"prob_pds_sir\n\nPositive and conservative autonomous nonlinear PDS\n\nbeginaligned\nu_1 = -2u_1u_2\nu_2 = 2u_1u_2 - u_2\nu_3 = u_2\nendaligned\n\nwith initial value mathbfu_0 = (099 0005 0005)^T and time domain (00 200). There is one independent linear invariant, e.g. u_1+u_2+u_3 = 10.\n\nReferences\n\nRonald E. Mickens, and Talitha M. Washington. \"NSFD discretizations of interacting population models satisfying conservation laws.\" Computers and Mathematics with Applications 66 (2013): 2307-2316. DOI: 10.1016/j.camwa.2013.06.011\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.prob_pds_stratreac","page":"API reference","title":"PositiveIntegrators.prob_pds_stratreac","text":"prob_pds_stratreac\n\nPositive and nonconservative autonomous nonlinear PDS\n\nbeginaligned\nu_1 = r_5 - r_6 -  r_7\nu_2 = 2r_1 - r_2 + r_3 - r_4 + r_6 - r_9 + r_10 - r_11\nu_3 = r_2 - r_3 - r_4 - r_5 - r_7 - r_8\nu_4 = -r_1 -r_2 + r_3 + 2r_4+r_5+2r_7+r_8+r_9\nu_5 = -r_8+r_9+r_10-r_11\nu_6 = r_8-r_9-r_10+r_11\nendaligned\n\nwith reaction rates\n\nbeginaligned\nr_1 =2643 10^-10σ^3 u_4  r_2 =801810^-17u_2 u_4   r_3 =61210^-4σ u_3\nr_4 =156710^-15u_3 u_2   r_5 = 107 10^-3σ^2u_3   r_6 = 71110^-11 81210^6 u_1\nr_7 = 1210^-10u_1 u_3  r_8 = 606210^-15u_3 u_5  r_9 = 106910^-11u_6 u_2\nr_10 = 128910^-2σ u_6  r_11 = 10^-8u_5 u_2\nendaligned\n\nwhere\n\nbeginaligned\nT = t3600 mod 24quad T_r=45quad T_s = 195\nσ(T) = begincases1  T_r T T_s0  textotherwiseendcases\nendaligned\n\nThe initial value is mathbfu_0 = (990610^1 662410^8 532610^11 169710^16 410^6 109310^9)^T and the time domain (432 10^4 302410^5). There are two independent linear invariants, e.g. u_1+u_2+3u_3+2u_4+u_5+2u_6=(113212)cdotmathbfu_0 and u_5+u_6 = 109710^9.\n\nReferences\n\nStephan Nüsslein, Hendrik Ranocha, and David I. Ketcheson. \"Positivity-preserving adaptive Runge-Kutta methods.\" Communications in Applied Mathematics and Computer Science 16 (2021): 155-179. DOI: 10.2140/camcos.2021.16.155\n\n\n\n\n\n","category":"constant"},{"location":"api_reference/#PositiveIntegrators.ConservativePDSProblem","page":"API reference","title":"PositiveIntegrators.ConservativePDSProblem","text":"ConservativePDSProblem(P, u0, tspan, p = NullParameters();\n                       p_prototype = nothing,\n                       analytic = nothing)\n\nA structure describing a conservative system of ordinary differential equation in form of a production-destruction system (PDS). P denotes the production matrix. u0 is the vector of initial conditions and tspan the time span (t_initial, t_final) of the problem. The optional argument p can be used to pass additional parameters to the function P.\n\nThe function P can be given either in the out-of-place form with signature production_terms = P(u, p, t) or the in-place form P(production_terms, u, p, t).\n\nKeyword arguments:\n\np_prototype: If P is given in in-place form, p_prototype is used to store evaluations of P.   If p_prototype is not specified explicitly and P is in-place, then p_prototype will be internally set to zeros(eltype(u0), (length(u0), length(u0))).\nanalytic: The analytic solution of a PDS must be given in the form f(u0,p,t).   Specifying the analytic solution can be useful for plotting and convergence tests.\n\nReferences\n\nHans Burchard, Eric Deleersnijder, and Andreas Meister. \"A high-order conservative Patankar-type discretisation for stiff systems of production-destruction equations.\" Applied Numerical Mathematics 47.1 (2003): 1-30. DOI: 10.1016/S0168-9274(03)00101-6\n\n\n\n\n\n","category":"type"},{"location":"api_reference/#PositiveIntegrators.MPE","page":"API reference","title":"PositiveIntegrators.MPE","text":"MPE([linsolve = ...])\n\nThe first-order modified Patankar-Euler algorithm for production-destruction systems. This one-step, one-stage method is first-order accurate, unconditionally positivity-preserving, and linearly implicit.\n\nThe scheme was introduced by Burchard et al for conservative production-destruction systems.  For nonconservative production–destruction systems we use the straight forward extension\n\nu_i^n+1 = u_i^n + Δt sum_j ji biggl(p_ij^n fracu_j^n+1u_j^n-d_ij^n fracu_i^n+1u_i^nbiggr) + Deltat p_ii^n - Δt d_ii^nfracu_i^n+1u_i^n.\n\nThe modified Patankar-Euler method requires the special structure of a PDSProblem or a ConservativePDSProblem.\n\nYou can optionally choose the linear solver to be used by passing an algorithm from LinearSolve.jl as keyword argument linsolve.\n\nReferences\n\nHans Burchard, Eric Deleersnijder, and Andreas Meister. \"A high-order conservative Patankar-type discretisation for stiff systems of production-destruction equations.\" Applied Numerical Mathematics 47.1 (2003): 1-30. DOI: 10.1016/S0168-9274(03)00101-6\n\n\n\n\n\n","category":"type"},{"location":"api_reference/#PositiveIntegrators.MPRK22","page":"API reference","title":"PositiveIntegrators.MPRK22","text":"MPRK22(α; [linsolve = ...])\n\nThe second-order modified Patankar-Runge-Kutta algorithm for  production-destruction systems. This one-step, two-stage method is second-order accurate, unconditionally positivity-preserving, and linearly implicit. The parameter α is described by Kopecz and Meister (2018) and studied by Izgin, Kopecz and Meister (2022) as well as Torlo, Öffner and Ranocha (2022).\n\nThe scheme was introduced by Kopecz and Meister for conservative production-destruction systems.  For nonconservative production–destruction systems we use the straight forward extension analogous to MPE.\n\nThis modified Patankar-Runge-Kutta method requires the special structure of a PDSProblem or a ConservativePDSProblem.\n\nYou can optionally choose the linear solver to be used by passing an algorithm from LinearSolve.jl as keyword argument linsolve.\n\nReferences\n\nHans Burchard, Eric Deleersnijder, and Andreas Meister. \"A high-order conservative Patankar-type discretisation for stiff systems of production-destruction equations.\" Applied Numerical Mathematics 47.1 (2003): 1-30. DOI: 10.1016/S0168-9274(03)00101-6\nStefan Kopecz and Andreas Meister. \"On order conditions for modified Patankar-Runge-Kutta schemes.\" Applied Numerical Mathematics 123 (2018): 159-179. DOI: 10.1016/j.apnum.2017.09.004\nThomas Izgin, Stefan Kopecz, and Andreas Meister. \"On Lyapunov stability of positive and conservative time integrators and application to second order modified Patankar-Runge-Kutta schemes.\" ESAIM: Mathematical Modelling and Numerical Analysis 56.3 (2022): 1053-1080. DOI: 10.1051/m2an/2022031\nDavide Torlo, Philipp Öffner, and Hendrik Ranocha. \"Issues with positivity-preserving Patankar-type schemes.\" Applied Numerical Mathematics 182 (2022): 117-147. DOI: 10.1016/j.apnum.2022.07.014\n\n\n\n\n\n","category":"type"},{"location":"api_reference/#PositiveIntegrators.MPRK43I","page":"API reference","title":"PositiveIntegrators.MPRK43I","text":"MPRK43I(α, β; [linsolve = ...])\n\nA family of third-order modified Patankar-Runge-Kutta schemes for (conservative) production-destruction systems, which is based on the two-parameter family of third order explicit Runge–Kutta schemes. Each member of this family is a one-step method with four-stages which is third-order accurate, unconditionally positivity-preserving, conservative and linearly implicit. In this implementation the stage-values are conservative as well. The parameters α and β must be chosen such that the Runge–Kutta coefficients are nonnegative,  see Kopecz and Meister (2018) for details. \n\nThe scheme was introduced by Kopecz and Meister for conservative production-destruction systems.  For nonconservative production–destruction systems we use the straight forward extension analogous to MPE.\n\nThese modified Patankar-Runge-Kutta methods require the special structure of a PDSProblem or a ConservativePDSProblem.\n\nYou can optionally choose the linear solver to be used by passing an algorithm from LinearSolve.jl as keyword argument linsolve.\n\nReferences\n\nStefan Kopecz and Andreas Meister. \"Unconditionally positive and conservative third order modified Patankar–Runge–Kutta   discretizations of production–destruction systems.\"  BIT Numerical Mathematics 58 (2018): 691–728. DOI: 10.1007/s10543-018-0705-1\n\n\n\n\n\n","category":"type"},{"location":"api_reference/#PositiveIntegrators.MPRK43II","page":"API reference","title":"PositiveIntegrators.MPRK43II","text":"MPRK43II(γ; [linsolve = ...])\n\nA family of third-order modified Patankar-Runge-Kutta schemes for (conservative) production-destruction systems, which is based on the one-parameter family of third order explicit Runge–Kutta schemes with  non-negative Runge–Kutta coefficients. Each member of this family is a one-step method with four stages which is third-order accurate, unconditionally positivity-preserving, conservative and linearly implicit. In this implementation the stage-values are conservative as well. The parameter γ must satisfy 3/8 ≤ γ ≤ 3/4.  Further details are given in Kopecz and Meister (2018).  \n\nThe scheme was introduced by Kopecz and Meister for conservative production-destruction systems.  For nonconservative production–destruction systems we use the straight forward extension analogous to MPE.\n\nThese modified Patankar-Runge-Kutta methods require the special structure of a PDSProblem or a ConservativePDSProblem.\n\nYou can optionally choose the linear solver to be used by passing an algorithm from LinearSolve.jl as keyword argument linsolve.\n\nReferences\n\nStefan Kopecz and Andreas Meister. \"Unconditionally positive and conservative third order modified Patankar–Runge–Kutta   discretizations of production–destruction systems.\"  BIT Numerical Mathematics 58 (2018): 691–728. DOI: 10.1007/s10543-018-0705-1\n\n\n\n\n\n","category":"type"},{"location":"api_reference/#PositiveIntegrators.PDSProblem","page":"API reference","title":"PositiveIntegrators.PDSProblem","text":"PDSProblem(P, D, u0, tspan, p = NullParameters();\n                   p_prototype = nothing,\n                   d_prototype = nothing,\n                   analytic = nothing)\n\nA structure describing a system of ordinary differential equations in form of a production-destruction system (PDS). P denotes the production matrix. The diagonal of P contains production terms without destruction counterparts. D is the vector of destruction terms without production counterparts. u0 is the vector of initial conditions and tspan the time span (t_initial, t_final) of the problem. The optional argument p can be used to pass additional parameters to the functions P and D.\n\nThe functions P and D can be used either in the out-of-place form with signature production_terms = P(u, p, t) or the in-place form P(production_terms, u, p, t).\n\nKeyword arguments:\n\np_prototype: If P is given in in-place form, p_prototype is used to store evaluations of P.   If p_prototype is not specified explicitly and P is in-place, then p_prototype will be internally set to zeros(eltype(u0), (length(u0), length(u0))).\nd_prototype: If D is given in in-place form, d_prototype is used to store evaluations of D. If d_prototype is not specified explicitly and D is in-place, then d_prototype will be internally\n\nset to zeros(eltype(u0), (length(u0),)).\n\nanalytic: The analytic solution of a PDS must be given in the form f(u0,p,t).   Specifying the analytic solution can be useful for plotting and convergence tests.\n\nReferences\n\nHans Burchard, Eric Deleersnijder, and Andreas Meister. \"A high-order conservative Patankar-type discretisation for stiff systems of production-destruction equations.\" Applied Numerical Mathematics 47.1 (2003): 1-30. DOI: 10.1016/S0168-9274(03)00101-6\n\n\n\n\n\n","category":"type"},{"location":"code_of_conduct/","page":"Code of conduct","title":"Code of conduct","text":"EditURL = \"https://github.com/SKopecz/PositiveIntegrators.jl/blob/main/CODE_OF_CONDUCT.md\"","category":"page"},{"location":"code_of_conduct/#code-of-conduct","page":"Code of conduct","title":"Code of Conduct","text":"","category":"section"},{"location":"code_of_conduct/","page":"Code of conduct","title":"Code of conduct","text":"Contributor Covenant Code of ConductOur PledgeWe as members, contributors, and leaders pledge to make participation in our community a harassment-free experience for everyone, regardless of age, body size, visible or invisible disability, ethnicity, sex characteristics, gender identity and expression, level of experience, education, socio-economic status, nationality, personal appearance, race, religion, or sexual identity and orientation.We pledge to act and interact in ways that contribute to an open, welcoming, diverse, inclusive, and healthy community.Our StandardsExamples of behavior that contributes to a positive environment for our community include:Demonstrating empathy and kindness toward other people\nBeing respectful of differing opinions, viewpoints, and experiences\nGiving and gracefully accepting constructive feedback\nAccepting responsibility and apologizing to those affected by our mistakes, and learning from the experience\nFocusing on what is best not just for us as individuals, but for the overall communityExamples of unacceptable behavior include:The use of sexualized language or imagery, and sexual attention or advances of any kind\nTrolling, insulting or derogatory comments, and personal or political attacks\nPublic or private harassment\nPublishing others' private information, such as a physical or email address, without their explicit permission\nOther conduct which could reasonably be considered inappropriate in a professional settingEnforcement ResponsibilitiesCommunity leaders are responsible for clarifying and enforcing our standards of acceptable behavior and will take appropriate and fair corrective action in response to any behavior that they deem inappropriate, threatening, offensive, or harmful.Community leaders have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct, and will communicate reasons for moderation decisions when appropriate.ScopeThis Code of Conduct applies within all community spaces, and also applies when an individual is officially representing the community in public spaces. Examples of representing our community include using an official e-mail address, posting via an official social media account, or acting as an appointed representative at an online or offline event.EnforcementInstances of abusive, harassing, or otherwise unacceptable behavior may be reported to Stefan Kopecz or Hendrik Ranocha. All complaints will be reviewed and investigated promptly and fairly.All community leaders are obligated to respect the privacy and security of the reporter of any incident.Enforcement GuidelinesCommunity leaders will follow these Community Impact Guidelines in determining the consequences for any action they deem in violation of this Code of Conduct:1. CorrectionCommunity Impact: Use of inappropriate language or other behavior deemed unprofessional or unwelcome in the community.Consequence: A private, written warning from community leaders, providing clarity around the nature of the violation and an explanation of why the behavior was inappropriate. A public apology may be requested.2. WarningCommunity Impact: A violation through a single incident or series of actions.Consequence: A warning with consequences for continued behavior. No interaction with the people involved, including unsolicited interaction with those enforcing the Code of Conduct, for a specified period of time. This includes avoiding interactions in community spaces as well as external channels like social media. Violating these terms may lead to a temporary or permanent ban.3. Temporary BanCommunity Impact: A serious violation of community standards, including sustained inappropriate behavior.Consequence: A temporary ban from any sort of interaction or public communication with the community for a specified period of time. No public or private interaction with the people involved, including unsolicited interaction with those enforcing the Code of Conduct, is allowed during this period. Violating these terms may lead to a permanent ban.4. Permanent BanCommunity Impact: Demonstrating a pattern of violation of community standards, including sustained inappropriate behavior,  harassment of an individual, or aggression toward or disparagement of classes of individuals.Consequence: A permanent ban from any sort of public interaction within the community.AttributionThis Code of Conduct is adapted from the [Contributor Covenant][homepage], version 2.0, available at https://www.contributor-covenant.org/version/2/0/codeofconduct.html.Community Impact Guidelines were inspired by Mozilla's code of conduct enforcement ladder.[homepage]: https://www.contributor-covenant.orgFor answers to common questions about this code of conduct, see the FAQ at https://www.contributor-covenant.org/faq. Translations are available at https://www.contributor-covenant.org/translations.","category":"page"},{"location":"license/","page":"License","title":"License","text":"EditURL = \"https://github.com/SKopecz/PositiveIntegrators.jl/blob/main/LICENSE.md\"","category":"page"},{"location":"license/#License","page":"License","title":"License","text":"","category":"section"},{"location":"license/","page":"License","title":"License","text":"MIT LicenseCopyright (c) 2023-present Stefan Kopecz, Hendrik Ranocha, and contributorsPermission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the \"Software\"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.","category":"page"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"EditURL = \"https://github.com/SKopecz/PositiveIntegrators.jl/blob/main/CONTRIBUTING.md\"","category":"page"},{"location":"contributing/#Contributing","page":"Contributing","title":"Contributing","text":"","category":"section"},{"location":"contributing/","page":"Contributing","title":"Contributing","text":"ContributingPositiveIntegrators.jl is an open-source project and we are very happy to accept contributions from the community. Please feel free to open issues or submit patches (preferably as pull requests) any time. For planned larger contributions, it is often beneficial to get in contact first, for example via issues.PositiveIntegrators.jl and its contributions are licensed under the MIT license (see License). As a contributor, you certify that all your contributions are in conformance with the Developer Certificate of Origin (Version 1.1), which is reproduced below.Developer Certificate of Origin (Version 1.1)The following text was taken from https://developercertificate.org:Developer Certificate of Origin\nVersion 1.1\n\nCopyright (C) 2004, 2006 The Linux Foundation and its contributors.\n1 Letterman Drive\nSuite D4700\nSan Francisco, CA, 94129\n\nEveryone is permitted to copy and distribute verbatim copies of this\nlicense document, but changing it is not allowed.\n\n\nDeveloper's Certificate of Origin 1.1\n\nBy making a contribution to this project, I certify that:\n\n(a) The contribution was created in whole or in part by me and I\n    have the right to submit it under the open source license\n    indicated in the file; or\n\n(b) The contribution is based upon previous work that, to the best\n    of my knowledge, is covered under an appropriate open source\n    license and I have the right under that license to submit that\n    work with modifications, whether created in whole or in part\n    by me, under the same open source license (unless I am\n    permitted to submit under a different license), as indicated\n    in the file; or\n\n(c) The contribution was provided directly to me by some other\n    person who certified (a), (b) or (c) and I have not modified\n    it.\n\n(d) I understand and agree that this project and the contribution\n    are public and that a record of the contribution (including all\n    personal information I submit with it, including my sign-off) is\n    maintained indefinitely and may be redistributed consistent with\n    this project or the open source license(s) involved.","category":"page"},{"location":"#PositiveIntegrators.jl","page":"Home","title":"PositiveIntegrators.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Julia library PositiveIntegrators.jl provides several time integration methods developed to preserve the positivity of numerical solutions.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"PositiveIntegrators.jl is a registered Julia package. Thus, you can install it from the Julia REPL via","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.add(\"PositiveIntegrators\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"If you want to update PositiveIntegrators.jl, you can use","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.update(\"PositiveIntegrators\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"As usual, if you want to update PositiveIntegrators.jl and all other packages in your current project, you can execute","category":"page"},{"location":"","page":"Home","title":"Home","text":"julia> using Pkg; Pkg.update()","category":"page"},{"location":"#Basic-examples","page":"Home","title":"Basic examples","text":"","category":"section"},{"location":"#Modified-Patankar-Runge-Kutta-schemes","page":"Home","title":"Modified Patankar-Runge-Kutta schemes","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Modified Patankar-Runge-Kutta (MPRK) schemes are unconditionally positive and conservative time integration schemes for the solution of positive and conservative ODE systems. The application of these methods is based on the representation of the ODE system as a so-called production-destruction system (PDS).","category":"page"},{"location":"#Production-destruction-systems-(PDS)","page":"Home","title":"Production-destruction systems (PDS)","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The application of MPRK schemes requires the ODE system to be represented as a production-destruction system (PDS). A PDS takes the general form","category":"page"},{"location":"","page":"Home","title":"Home","text":"    u_i(t) = sum_j=1^N bigl(p_ij(tboldsymbol u) - d_ij(tboldsymbol u)bigr)quad i=1dotsN","category":"page"},{"location":"","page":"Home","title":"Home","text":"where boldsymbol u=(u_1dotsu_n)^T is the vector of unknowns and both production terms p_ij(tboldsymbol u) and destruction terms d_ij(tboldsymbol u) must be nonnegative for all ij=1dotsN. The meaning behind p_ij and d_ij is as follows:","category":"page"},{"location":"","page":"Home","title":"Home","text":"p_ij with ine j represents the sum of all nonnegative terms which appear in equation i with a positive sign and in equation j with a negative sign.\nd_ij with ine j represents the sum of all nonnegative terms which appear in equation i with a negative sign and in equation j with a positive sign.\np_ii represents the sum of all nonnegative terms  which appear in equation i and don't have a negative counterpart in one of the other equations.\nd_ii represents the sum of all negative terms which appear in equation i and don't have a positive counterpart in one of the other equations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"This naming convention leads to p_ij = d_ji for i j and therefore a PDS is completely defined by the production matrix mathbfP=(p_ij)_ij=1dotsN and the destruction vector mathbfd=(d_ii)_i=1dotsN.","category":"page"},{"location":"","page":"Home","title":"Home","text":"As an example we consider the Lotka-Volterra model","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\nu_1 = 2u_1-u_1u_2\nu_2 = u_1u_2-u_2\nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"which always has positive solutions if positive initial values are supplied. Assuming u_1u_20, the above naming scheme results in","category":"page"},{"location":"","page":"Home","title":"Home","text":"beginaligned\np_11(u_1u_2) = 2u_1\np_21(u_1u_2) = u_1u_2 = d_12(u_1u_2) \nd_22(u_1u_2) = u_2\nendaligned","category":"page"},{"location":"","page":"Home","title":"Home","text":"where all remaining production and destruction terms are zero. Consequently the production matrix mathbf P and destruction vector mathbf d are","category":"page"},{"location":"","page":"Home","title":"Home","text":"mathbf P(u_1u_2) = beginpmatrix2u_1  0 u_1u_2  0endpmatrixquad mathbf d(u_1u_2) = beginpmatrix0 u_2endpmatrix","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg; Pkg.add(\"OrdinaryDiffEq\");  Pkg.add(\"Plots\")","category":"page"},{"location":"","page":"Home","title":"Home","text":"To solve this PDS together with initial values u_1(0)=u_2(0)=2 on the time domain (010), we first need to create a PDSProblem.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using PositiveIntegrators # load PDSProblem\n\nP(u, p, t) = [2*u[1]  0.0; u[1]*u[2]  0.0] # Production matrix\nd(u, p, t) = [0.0; u[2]] # Destruction vector\n\nu0 = [2.0; 2.0] # initial values\ntspan = (0.0, 10.0) # time span\n\n# Create PDS\nprob = PDSProblem(P, d, u0, tspan)\nnothing #hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"Now that the problem has been created, we can solve it with any of the methods of OrdinaryDiffEq.jl. Here we use the method Tsit5(). Please note that PositiveIntegrators.jl currently only provides methods for positive and conservative PDS, see below.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using OrdinaryDiffEq  #load Tsit5\n\nsol = solve(prob, Tsit5())\nnothing # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, we can use Plots.jl to visualize the solution.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Plots\n\nplot(sol)","category":"page"},{"location":"#Conservative-production-destruction-systems","page":"Home","title":"Conservative production-destruction systems","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A PDS with the additional property","category":"page"},{"location":"","page":"Home","title":"Home","text":"  p_ii(tboldsymbol y)=d_ii(tboldsymbol y)=0","category":"page"},{"location":"","page":"Home","title":"Home","text":"for i=1dotsN is called conservative. In this case we have p_ij=d_ji for all ij=1dotsN, which leads to","category":"page"},{"location":"","page":"Home","title":"Home","text":"fracddtsum_i=1^N y_i=sum_i=1^N y_i = sum_mathclapij=1^N bigl(p_ij(tboldsymbol y) - d_ij(tboldsymbol y)bigr)= sum_mathclapij=1^N bigl(p_ij(tboldsymbol y) - p_ji(tboldsymbol y)bigr) = 0","category":"page"},{"location":"","page":"Home","title":"Home","text":"This shows that the sum of the state variables of a conservative PDS remains constant over time, i.e.","category":"page"},{"location":"","page":"Home","title":"Home","text":"sum_i=1^N y_i(t) = sum_i=1^N y_i(0)","category":"page"},{"location":"","page":"Home","title":"Home","text":"for all times t0. Moreover, a conservative PDS is completely defined by the square matrix mathbf P=(p_ij)_ij=1dotsN. There is no need to store an additional vector of destruction terms since d_ij = p_ji for all ij=1dotsN.","category":"page"},{"location":"","page":"Home","title":"Home","text":"One specific example of a conservative PDS is the SIR model","category":"page"},{"location":"","page":"Home","title":"Home","text":"S = -fracβ S INquad I= fracβ S IN - γ Iquad R=γ I","category":"page"},{"location":"","page":"Home","title":"Home","text":"with N=S+I+R and betagamma0. Assuming SIR0 the production and destruction terms are given by","category":"page"},{"location":"","page":"Home","title":"Home","text":"p_21(SIR) = d_12(SIR) = fracβ S INquad p_32(SIR) = d_23(SIR) = γ I","category":"page"},{"location":"","page":"Home","title":"Home","text":"where the remaining production and destruction terms are zero. The corresponding production matrix mathbf P is","category":"page"},{"location":"","page":"Home","title":"Home","text":"mathbf P(SIR) = beginpmatrix0  0  0 fracβ S IN  0  0 0  γ I  0endpmatrix","category":"page"},{"location":"","page":"Home","title":"Home","text":"The following example shows how to implement the above SIR model with beta=04 gamma=004, initial conditions S(0)=997 I(0)=3 R(0)=0 and time domain (0 100) using ConservativePDSProblem from PositiveIntegrators.jl.","category":"page"},{"location":"","page":"Home","title":"Home","text":"import Pkg; Pkg.add(\"OrdinaryDiffEq\");","category":"page"},{"location":"","page":"Home","title":"Home","text":"using PositiveIntegrators\n\n# Out-of-place implementation of the P matrix for the SIR model\nfunction P(u, p, t)\n  S, I, R = u\n\n  β = 0.4\n  γ = 0.04\n  N = 1000.0\n\n  P = zeros(3,3)\n  P[2,1] = β*S*I/N\n  P[3,2] = γ*I\n  return P\nend\n\nu0 = [997.0; 3.0; 0.0]; # initial values\ntspan = (0.0, 100.0); # time span\n\n# Create SIR problem\nprob = ConservativePDSProblem(P, u0, tspan)\nnothing # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"Since the SIR model is not only conservative but also positive, we can use any MPRK scheme from PositiveIntegrators.jl to solve it. Here we use MPRK22(1.0). Please note that any method from OrdinaryDiffEq.jl can be used as well, but might possibly generate negative approximations.","category":"page"},{"location":"","page":"Home","title":"Home","text":"sol = solve(prob, MPRK22(1.0))\nnothing # hide","category":"page"},{"location":"","page":"Home","title":"Home","text":"Finally, we can use Plots.jl to visualize the solution.","category":"page"},{"location":"","page":"Home","title":"Home","text":"using Plots\n\nplot(sol, legend=:right)","category":"page"},{"location":"#Referencing","page":"Home","title":"Referencing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"If you use PositiveIntegrators.jl for your research, please cite it using the bibtex entry","category":"page"},{"location":"","page":"Home","title":"Home","text":"@misc{PositiveIntegrators.jl,\n  title={{PositiveIntegrators.jl}: {A} {J}ulia library of positivity-preserving\n         time integration methods},\n  author={Kopecz, Stefan and Ranocha, Hendrik and contributors},\n  year={2023},\n  doi={10.5281/zenodo.10868393},\n  url={https://github.com/SKopecz/PositiveIntegrators.jl}\n}","category":"page"},{"location":"#License-and-contributing","page":"Home","title":"License and contributing","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"This project is licensed under the MIT license (see License). Since it is an open-source project, we are very happy to accept contributions from the community. Please refer to the section Contributing for more details.","category":"page"}]
}

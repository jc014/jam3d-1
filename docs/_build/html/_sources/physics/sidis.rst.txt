SIDIS
=====

The differential cross section is given by 

.. math::

  \begin{align}
  \frac{d\sigma}{dx ~ dy ~ d\phi_S ~ dz ~ d\phi_h ~ dP_{\rm T}^2} = 
  \frac{\alpha^2}{xyQ^2}\frac{y^2}{2(1-\varepsilon)}\left(1+\frac{\gamma^2}{2x} \right) 
  \sum_{i=1}^{18} F_i(x,z,Q^2,P_{\rm T}) \beta_i 
  \end{align}

where 

.. list-table:: SIDIS structure functions
   :widths: 5 5 5
   :header-rows: 1

   * - :math:`F_i`
     - Standard Label
     - :math:`\beta_i`
   * - :math:`F_1`
     - :math:`F_{UU,T}`
     - :math:`1`
   * - :math:`F_2`
     - :math:`F_{UU,L}`
     - :math:`\varepsilon`
   * - :math:`F_3`
     - :math:`F_{LL}`
     - :math:`S_{||}\lambda_e\sqrt{1-\varepsilon^2}`
   * - :math:`F_4`
     - :math:`F_{UT}^{\sin(\phi_h+\phi_S)}`
     - :math:`|\vec{S}_\perp|\varepsilon\sin(\phi_h+\phi_S)`
   * - :math:`F_5`
     - :math:`F_{UT,T}^{\sin(\phi_h-\phi_S)}`
     - :math:`|\vec{S}_\perp|\sin(\phi_h-\phi_S)`
   * - :math:`F_6`
     - :math:`F_{UT,L}^{\sin(\phi_h-\phi_S)}`
     - :math:`|\vec{S}_\perp|\varepsilon\sin(\phi_h-\phi_S)`
   * - :math:`F_7`
     - :math:`F_{UU}^{\cos2\phi_h}`
     - :math:`\varepsilon~\cos(2\phi_h)`
   * - :math:`F_8`
     - :math:`F_{UT}^{\sin(3\phi_h-\phi_S)}`
     - :math:`|\vec{S}_\perp|\varepsilon\sin(3\phi_h-\phi_S)`
   * - :math:`F_9`
     - :math:`F_{LT}^{\cos(\phi_h-\phi_S)}`
     - :math:`|\vec{S}_\perp|\lambda_e\sqrt{1-\varepsilon^2}\cos(\phi_h-\phi_S)`
   * - :math:`F_{10}`
     - :math:`F_{UL}^{\sin~2\phi_h}`
     - :math:`S_{||}\varepsilon\sin(2\phi_h)`
   * - :math:`F_{11}`
     - :math:`F_{LT}^{\cos\phi_S}`
     - :math:`|\vec{S}_\perp|\lambda_e\sqrt{2\varepsilon(1-\varepsilon)}\cos\phi_S`
   * - :math:`F_{12}`
     - :math:`F_{LL}^{\cos\phi_h}`
     - :math:`S_{||}\lambda_e\sqrt{2\varepsilon(1-\varepsilon)}\cos\phi_h`
   * - :math:`F_{13}`
     - :math:`F_{LT}^{\cos(2\phi_h-\phi_S)}`
     - :math:`|\vec{S}_\perp|\lambda_e\sqrt{2\varepsilon(1-\varepsilon)}\cos(2\phi_h-\phi_S)`
   * - :math:`F_{14}`
     - :math:`F_{UL}^{\sin\phi_h}`
     - :math:`S_{||}\sqrt{2\varepsilon(1+\varepsilon)}\sin\phi_h`
   * - :math:`F_{15}`
     - :math:`F_{LU}^{\sin\phi_h}`
     - :math:`\lambda_e\sqrt{2\varepsilon(1-\varepsilon)}\sin\phi_h`
   * - :math:`F_{16}`
     - :math:`F_{UU}^{\cos\phi_h}`
     - :math:`\sqrt{2\varepsilon(1+\varepsilon)}~\cos\phi_h`
   * - :math:`F_{17}`
     - :math:`F_{UT}^{\sin\phi_S}`
     - :math:`|\vec{S}_\perp|\sqrt{2\varepsilon(1+\varepsilon)}\sin\phi_S`
   * - :math:`F_{18}`
     - :math:`F_{UT}^{\sin(2\phi_h-\phi_S)}`
     - :math:`|\vec{S}_\perp|\sqrt{2\varepsilon(1+\varepsilon)}\sin(2\phi_h-\phi_S)`


were 

.. math::

  \begin{align}
  \varepsilon &= \frac{1-y-\frac{1}{4}\gamma^2 y^2}{1-y+\frac{1}{2} y^2 + \frac{1}{4}\gamma^2 y^2}\, , \\
  \gamma &=  \frac{2 M x}{Q} \, ,
  \end{align}

and in the nucleon rest frame the polarization vector is given by :math:`S=(0,\vec{S}_\perp,S_L)` with :math:`\vec{S}_\perp^2+S_L^2=1`.

The 18 structure function in SIDIS at leading-order will be expressed
in the context of WW-type approximation in terms of a minimal
TMD basis using the gaussian ansatz:

.. math::

  \begin{align}
  {\cal F}_q(\xi,p_{\perp})&={\cal K}_q ~ {\cal C}_q(\xi) \frac{\exp\left(-k_{\perp}^2/\omega_q\right)}{\pi \omega_q}\\
  {\cal D}_q(\xi,p_{\perp})&={\cal K}_q ~ {\cal C}_q(\xi) \frac{\exp\left(-P_{\perp}^2/\omega_q\right)}{\pi \omega_q}.
  \end{align}

We denote the transverse momentum of the quark inside a fast moving
proton by :math:`k`. We use the notation :math:`p_{\perp}` for the
transverse momentum of the quark relative to the fragmenting parton
motion. The structure functions are expressed as

.. math::

  \begin{align}
  F&=\sum_q e_q^2 ~ {\cal K}_q ~ {\cal F}_q(x) ~{\cal D}_q(z) \frac{\exp\left(-P_{\perp}^2/\Omega_q\right)}{\pi \Omega_q}\\
  \Omega_q&=z^2\left<k\right> + \left<p_{\perp}\right>
  \end{align}

.. list-table:: Nucleon distributions and fragmentation functions
   :header-rows: 1

   * - type
     - Name
     - :math:`{\cal K}_q`
     - :math:`{\cal C}_q`
   * - :math:`{\cal F}_q`
     - upol.PDF
     - :math:`1`
     - :math:`f_1`
   * - :math:`{\cal F}_q`
     - pol.PDF
     - :math:`1`
     - :math:`g_1^q`
   * - :math:`{\cal F}_q`
     - Transversity
     - :math:`1`
     - :math:`{h_1^q}`
   * - :math:`{\cal F}_q`
     - Sivers
     - :math:`\frac{2M^2}{\omega_q}`
     - :math:`{{f^{\perp(1)q}_{1T}}}`
   * - :math:`{\cal F}_q`
     - Boer-Mulders
     - :math:`\frac{2M^2}{\omega_q}`
     - :math:`{{h^{\perp(1)q}_1 }}`
   * - :math:`{\cal F}_q`
     - Pretzelosity
     - :math:`\frac{2M^4}{\omega_q^2}`
     - :math:`{{h^{\perp(2)q}_{1T} }}`
   * - :math:`{\cal F}_q`
     - WormGear
     - :math:`\frac{2M^2}{\omega_q}`
     - :math:`{{g^{\perp (1) q}_{1T}}}`
   * - :math:`{\cal F}_q`
     - WormGear
     - :math:`\frac{2M^2}{\omega_q}`
     - :math:`{{h^{\perp (1) q}_{1L}}}`
   * - :math:`{\cal F}_q`
     - twist-3
     - :math:`1`
     - :math:`{{g^{q}_{T}}}`
   * - :math:`{\cal F}_q`
     - twist-3
     - :math:`\frac{2M^4}{\omega_q^2}`
     - :math:`{ {g^{\perp(2)q}_{T} }}`
   * - :math:`{\cal F}_q`
     - twist-3
     - :math:`\frac{2M^2}{\omega_q}`
     - :math:`{{f^{\perp (1) q}}}`
   * - :math:`{\cal F}_q`
     - twist-3
     - :math:`\frac{2M^2}{\omega_q}`
     - :math:`{{h^{\perp (1) q}_T}}`
   * - :math:`{\cal F}_q`
     - twist-3
     - :math:`\frac{2M^2}{\omega_q}`
     - :math:`{{h^{(1) q}_T}}`
   * - :math:`{\cal F}_q`
     - twist-3
     - :math:`\frac{2M^4}{\omega_q^2}`
     - :math:`{{f^{\perp(2)q}_{T} }}`
   * - :math:`{\cal C}_q`
     - FF
     - :math:`1`
     - :math:`{{d_1^{q}}}`
   * - :math:`{\cal C}_q`
     - Collins
     - :math:`\frac{2z^2m_h^2}{\omega_q}`
     - :math:`{{H^{\perp(1)q}_{1} }}`


twist-3 functions can be related to twist-2 functions using the following formulas:

.. math::
   \begin{align}
   g^{q}_{T}(x) = \int_x^1\frac{ d y}{y}\,g_1^q(y) \\
   \end{align}

.. math::
   \begin{align}
   xg^{\perp(2)q}_{T}(x)  = \frac{\left<k_{\perp}^2\right>_{g^{\perp (1) q}_{1T}}}{M^2}\; g^{\perp (1) q}_{1T}(x)\\
   \end{align}


.. math::
   \begin{align}
   x\,f^{\perp (1) q}(x) = \frac{\left<k_{\perp}^2\right>_{f_1}}{2M^2}\,f_1(x)\\
   \end{align}

.. math::
   \begin{align}
   -\frac{1}{2}x(h^{\perp (1) q}_T(x)-h^{(1) q}_T(x)) = \frac{\left<k_{\perp}^2\right>_{h_1^q}}{2M^2}\;h_1^q(x)
   \end{align}

.. math::
   \begin{equation}
    -\frac{1}{2}x(h^{\perp (1) q}_T(x)+h^{(1) q}_T(x)) = h^{\perp(2)q}_{1T}(x)
   \end{equation}


.. list-table:: Structure functions in terms of partonic d.o.f
   :header-rows: 1

   * - :math:`F_i`
     - standard notation 
     - :math:`{\cal K}_q`
     - :math:`{\cal F}_q(x)`
     - :math:`{\cal D}_q(z)`
   * - :math:`F_1`
     - :math:`F_{UU,T}`
     - :math:`x`
     - :math:`f_1`
     - :math:`d_1`
   * - :math:`F_2`
     - :math:`F_{UU,L}`
     - :math:`0`
     - :math:`0`
     - :math:`0`
   * - :math:`F_3`
     - :math:`F_{LL}`
     - :math:`x`
     - :math:`{g_1^q}`
     - :math:`{d_1^q}`
   * - :math:`F_4`
     - :math:`F_{UT}^{\sin(\phi_h+\phi_S)}`
     - :math:`\frac{2xz{P^h_{\perp}} m_h}{w_q}`
     - :math:`{h_1^q}`
     - :math:`{H^{\perp(1)q}_{1}}`
   * - :math:`F_5`
     - :math:`F_{UT,T}^{\sin(\phi_h-\phi_S)}`
     - :math:`-\frac{2xzM{P^h_{\perp}}}{w_q}`
     - :math:`{f^{\perp(1)q}_{1T}}`
     - :math:`{d_1^q}`
   * - :math:`F_6`
     - :math:`F_{UT,L}^{\sin(\phi_h-\phi_S)}`
     - :math:`0`
     - :math:`0`
     - :math:`0`
   * - :math:`F_7`
     - :math:`F_{UU}^{\cos(2\phi_h)}`
     - :math:`\frac{4xz^2M{P^h_{\perp}}^2m_h}{w_q^2}`
     - :math:`{h^{\perp(1)q}_1}`
     - :math:`{H^{\perp(1)q}_{1}}`
   * - :math:`F_8`
     - :math:`F_{UT}^{\sin(3\phi_h-\phi_S)}`
     - :math:`\frac{2xz^3{P^h_{\perp}}^3m_hM^2}{w_q^3}`
     - :math:`{h^{\perp(2)q}_{1T}}`
     - :math:`{H^{\perp(1)q}_{1}}`
   * - :math:`F_9`
     - :math:`F_{LT}^{\cos(\phi_h-\phi_S)}`
     - :math:`\frac{2xzM{P^h_{\perp}}}{w_q}`
     - :math:`{g^{\perp (1) q}_{1T}}`
     - :math:`{d_1^q}`
   * - :math:`F_{10}`
     - :math:`F_{UL}^{\sin(2\phi_h)}`
     - :math:`\frac{4xz^2M{P^h_{\perp}}^2m_h}{w_q^2}`
     - :math:`{h^{\perp (1) q}_{1L}}`
     - :math:`{H^{\perp(1)q}_{1}}`
   * - :math:`F_{11}`
     - :math:`F_{LT}^{\cos\phi_S}`
     - :math:`-\frac{2M}{Q}x`
     - :math:`{g^{q}_{T}}`
     - :math:`{d_1^q}`
   * - :math:`F_{12}`
     - :math:`F_{LL}^{\cos\phi_h}`
     - :math:`-\frac{2xz{P^h_{\perp}}}{Q}\frac{{\left<p_{\perp}\right>}}{w_q}`
     - :math:`{g_1^q}`
     - :math:`{d_1^q}`
   * - :math:`F_{13}`
     - :math:`F_{LT}^{\cos(2\phi_h-\phi_S)}`
     - :math:`-\frac{2xz^2M^3{P^h_{\perp}}^2}{Q}\frac{1}{w_q^2}`
     - :math:`x{g^{q}_{T}}_{\perp}`
     - :math:`{d_1^q}`
   * - :math:`F_{14}`
     - :math:`F_{UL}^{\sin\phi_h}`
     - :math:`\frac{4Mm_h}{Q}\frac{{P^h_{\perp}}z}{w_q}`
     - :math:`{h^{\perp (1) q}_{1L}}`
     - :math:`{H^{\perp(1)q}_{1}}`
   * - :math:`F_{15}`
     - :math:`F_{LU}^{\sin\phi_h}`
     - :math:`0`
     - :math:`0`
     - :math:`0`
   * - :math:`F_{16}`
     - :math:`F_{UU}^{\cos\phi_h}(i)`
     - :math:`0`
     - :math:`{h^{\perp(1)q}_1 }`
     - :math:`{H^{\perp(1)q}_{1}}`
   * - :math:`F_{16}`
     - :math:`F_{UU}^{\cos\phi_h}(ii)`
     - :math:`-\frac{2M}{Q}\frac{2xz{P^h_{\perp}}M}{w_q}`
     - :math:`xf^{\perp (1) q}`
     - :math:`{d_1^q}`
   * - :math:`F_{17}`
     - :math:`F_{UT}^{\sin\phi_S}(i)`
     - :math:`0`
     - :math:`{f^{\perp(1)q}_{1T}}`
     - :math:`{d_1^q}`
   * - :math:`F_{17}`
     - :math:`F_{UT}^{\sin\phi_S}(ii)`
     - :math:`\frac{2M}{Q}\frac{4xz^2m_hM}{w_q}(1-\frac{{P^h_{\perp}}^2}{w_q})`
     - :math:`-\frac{1}{2}x({h^{\perp (1) q}_T}-{f^{\perp(2)q}_{T} })`
     - :math:`{H^{\perp(1)q}_{1}}`
   * - :math:`F_{18}`
     - :math:`F_{UT}^{\sin(2\phi_h-\phi_S)}(i)`
     - :math:`\frac{2M}{Q}x\frac{M^2z^2{P^h_{\perp}}^2}{w_q^2}`
     - :math:`x{f^{\perp (1) q}}`
     - :math:`{d_1^q}`
   * - :math:`F_{18}`
     - :math:`F_{UT}^{\sin(2\phi_h-\phi_S)}(ii)`
     - :math:`-\frac{2M^2}{Q}x\frac{4z^2{P^h_{\perp}}^2Mm_h}{w_q^2}`
     - :math:`-\frac{1}{2}x({h^{\perp (1) q}_T}+{f^{\perp(2)q}_{T} })`
     - :math:`{H^{\perp(1)q}_{1}}`
   * - :math:`F_{19}`
     - :math:`F_{\rm CAHN}^{\cos(2\phi_h)}`
     - :math:`\frac{1}{Q^2}\frac{2xz^2{P^h_{\perp}}^2{\left<p_{\perp}\right>}^2}{w_q^2}`
     - :math:`{f_1^q}`
     - :math:`{d_1^q}`







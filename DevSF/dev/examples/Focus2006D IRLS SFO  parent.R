guitest <- mkinmod.full(
  parent = list(
    sink  = TRUE,
    type = "SFO",
    to='m1',
    k = list(ini   = 0.040,
             fixed = 0,
             lower = 0.0,
             upper = Inf),
    M0 = list(ini   = 100.15,
              fixed = 0,
              lower = 0.0,
              upper = Inf)),
  m1=list(type = "SFO"),
  inpartri='default',outpartri='default',data=mkin_long_to_wide(FOCUS_2006_D) )

#
# Fit and optimizer
#

Fit    <- IRLSkinfit.full(
  guitest,
  plot      = TRUE,
  quiet     = TRUE,
  ctr       = kingui.control(
    method = 'solnp',
    submethod = 'Port',
    maxIter = 100,
    tolerance = 1E-06,
    odesolver = 'lsoda'),
  irls.control = list(
    maxIter = 10,
    tolerance = 0.001))

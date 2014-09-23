SET ( CF_KERNEL_MODS logcpp Common Config Environment MathTools Framework ShapeFunctions Solver )

FOREACH( ADIR ${CF_KERNEL_MODS} )
    LOGVERBOSE ("\# KERNEL MODULE [${ADIR}]")
    ADD_SUBDIRECTORY (${ADIR})
ENDFOREACH( ADIR ${CF_KERNEL_MODS} )

ADD_SUBDIRECTORY ( UnitTests )
LOGVERBOSE ( "\#  ADDING UNITTEST MODULE [UnitTests]")

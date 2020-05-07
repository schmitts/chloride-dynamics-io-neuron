TITLE Simplified NpHR channel for activation of Chloride currents at set levels

COMMENT
-----------------------------------------------------------------------------
	The inward chloride current is set as a persistent value for the duration
	of the stimulation.

	Authors: Christopher Brian Currin

-----------------------------------------------------------------------------
ENDCOMMENT

NEURON {
	THREADSAFE
	POINT_PROCESS NpHR
	USEION cl WRITE icl VALENCE -1
	RANGE del, dur, amp
	:ELECTRODE_CURRENT icl
}

UNITS {
	(nA) = (nanoamp)
}

PARAMETER {
	del = 0 (ms)
	dur = 0 (ms)
	amp = 0	(nA)
}

ASSIGNED {
	icl		(nA)	: chloride current
	on		(1)
}

INITIAL { 
	icl = 0
	on = 0
  	net_send(del, 1)
}

BREAKPOINT {
	if(on){
		icl  = amp
	}
	else{
		icl = 0
	}
}

NET_RECEIVE (w) {
  if (flag == 1) {
    if (on == 0) {
      : turn it on
      on = 1
      : prepare to turn it off
      net_send(dur, 1)
    } else {
      : turn it off
      on = 0
    }
  }
}

COMMENT

http://www.neuron.yale.edu/phpbb/viewtopic.php?f=16&t=2889

ENDCOMMENT
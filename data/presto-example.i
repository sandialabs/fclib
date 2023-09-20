## Copyright (2000) Sandia Corporation. Under the terms of Contract
## DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains 
## certain rights in this software.
##
## Redistribution and use in source and binary forms, with or without
## modification, are permitted provided that the following conditions are met:
##
##   * Redistributions of source code must retain the above copyright
##     notice, this list of conditions and the following disclaimer.
##
##   * Redistributions in binary form must reproduce the above copyright
##     notice, this list of conditions and the following disclaimer in the
##     documentation and/or other materials provided with the
##     distribution.
##
##   * Neither the name of Sandia nor the names of any contributors may
##     be used to endorse or promote products derived from this software
##     without specific prior written permission.
##
##   * Modified source versions must be plainly marked as such, and must
##     not be misrepresented as being the original software.
##
## THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
## AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
## IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
## ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE FOR 
## ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
## DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
## SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
## CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT 
## LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
## OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH 
## DAMAGE.

## $Source$
## $Revision$
## $Date$
#
# Generated by:
# SIMBA version 22 Build number 1164(Release).
# Built on sid at Fri Apr 4 09:27:11 PST 2003
#
begin sierra simple-assy-new-bolts
  
  begin property specification for material putty-plh-mm-g-s-2
    density = 7.8e-9
    begin parameters for model ep_power_hard
      hardening constant = 70
      hardening exponent = .1
      luders strain = 0
      poissons ratio = .3
      yield stress = 70
      youngs modulus = 210e3
    end parameters for model ep_power_hard
  end property specification for material putty-plh-mm-g-s-2
  
  begin property specification for material steel-elastic-mm-g-s
    density = 7.8e-9
    begin parameters for model elastic
      poissons ratio = .3
      youngs modulus = 210e3
    end parameters for model elastic
  end property specification for material steel-elastic-mm-g-s
  
  # Functions for conditions
  
  # Function spot-weld-normal
  begin definition for function spot-weld-normal
    type = piecewise linear
    begin values
      0. 0.
      .25 500.
      .5 1000.
      1. 2000.
      1.1 0.
    end values
  end definition for function spot-weld-normal
  
  # Function spot-weld-normal-2
  begin definition for function spot-weld-normal-2
    type = piecewise linear
    begin values
      0. 0.
      .25 1500.
      .5 3000.
      1. 6000.
      1.1 0.
    end values
  end definition for function spot-weld-normal-2
  
  # Function spot-weld-tangential
  begin definition for function spot-weld-tangential
    type = piecewise linear
    begin values
      0. 0.
      .25 250.
      .5 500.
      1. 1000.
      1.1 0.
    end values
  end definition for function spot-weld-tangential
  
  # Function spot-weld-tangential-2
  begin definition for function spot-weld-tangential-2
    type = piecewise linear
    begin values
      0. 0.
      .25 750.
      .5 1500.
      1. 3000.
      1.1 0.
    end values
  end definition for function spot-weld-tangential-2
  
  begin finite element model simple-assy-new-bolts
    database name = simple-assy-new-bolts.g
    database type = exodusII
    
    # - Block id 11 had name 11
    begin parameters for block block_11
      material steel-elastic-mm-g-s
      solid mechanics use model elastic
      hex hourglass stiffness = 0.05
      hex hourglass viscosity = 0.03
    end parameters for block block_11
    
    # - Block id 21 had name 21
    begin parameters for block block_21
      material steel-elastic-mm-g-s
      solid mechanics use model elastic
      hex hourglass stiffness = 0.05
      hex hourglass viscosity = 0.03
    end parameters for block block_21
    
    # - Block id 3 had name case
    begin parameters for block block_3
      material putty-plh-mm-g-s-2
      solid mechanics use model ep_power_hard
      hex hourglass stiffness = 0.05
      hex hourglass viscosity = 0.03
    end parameters for block block_3
    
    # - Block id 1 had name cover
    begin parameters for block block_1
      material putty-plh-mm-g-s-2
      solid mechanics use model ep_power_hard
      hex hourglass stiffness = 0.05
      hex hourglass viscosity = 0.03
    end parameters for block block_1
    
    # - Block id 4 had name floor
    begin parameters for block block_4
      material steel-elastic-mm-g-s
      solid mechanics use model elastic
    end parameters for block block_4
  end finite element model simple-assy-new-bolts
  
  # Direction is used by more than one block.
  define direction simba_v_shared_1 with vector 0. 0. -10000
  
  begin presto procedure Apst_Procedure
    
    #
    # *** Time step control information
    begin time control
      termination time = 1.e-3
      
      begin time stepping block time_stepping_block
        start time = 0
        begin parameters for presto region presto
          step interval = 100
        end parameters for presto region
      end time stepping block
    end time control
    
    begin presto region presto
      
      #
      # **** Results section.
      begin results output results
        database name = simple-assy-new-bolts.e
        database type = exodusII
        at time 0., increment = 50e-6
        title - test model
        nodal variables = displacement as displ
        nodal variables = velocity as vel
        nodal variables = acceleration as accel
        nodal variables = reactions as react
        nodal variables = force_external as f_ext
        nodal variables = force_internal as f_int
        nodal variables = force_contact as f_contact
        
        element variables = stress
        element variables = state_ep_power_hard as state_eph
        
        global variables = timestep
        global variables = kineticenergy as KE
      end results output results
      
      use finite element model simple-assy-new-bolts
      #
      # *** Contact section.
      begin contact definition contacts
        begin defaults
          friction model = sliding-no-friction
          normal tolerance = .005
          overlap normal tolerance = .005
          overlap tangential tolerance = .005
          tangential tolerance = .005
        end defaults
        
        #
        begin frictionless model sliding-no-friction
        end
        
        contact all blocks
        remove initial overlap
      end contact definition contacts
      #
      # *** Spot welds.
      #
      # Spot weld sw1
      begin spot weld sw1
        node set = nodelist_1010
        surface = surface_1000
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal-2
        tangential displacement function = spot-weld-tangential-2
      end spot weld sw1
      #
      # Spot weld sw10
      begin spot weld sw10
        node set = nodelist_1013
        surface = surface_1009
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal-2
        tangential displacement function = spot-weld-tangential-2
      end spot weld sw10
      #
      # Spot weld sw2
      begin spot weld sw2
        node set = nodelist_1014
        surface = surface_1001
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal
        tangential displacement function = spot-weld-tangential
      end spot weld sw2
      #
      # Spot weld sw3
      begin spot weld sw3
        node set = nodelist_1015
        surface = surface_1002
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal
        tangential displacement function = spot-weld-tangential
      end spot weld sw3
      #
      # Spot weld sw4
      begin spot weld sw4
        node set = nodelist_1016
        surface = surface_1003
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal
        tangential displacement function = spot-weld-tangential
      end spot weld sw4
      #
      # Spot weld sw5
      begin spot weld sw5
        node set = nodelist_1011
        surface = surface_1004
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal-2
        tangential displacement function = spot-weld-tangential-2
      end spot weld sw5
      #
      # Spot weld sw6
      begin spot weld sw6
        node set = nodelist_1012
        surface = surface_1005
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal-2
        tangential displacement function = spot-weld-tangential-2
      end spot weld sw6
      #
      # Spot weld sw7
      begin spot weld sw7
        node set = nodelist_1017
        surface = surface_1006
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal
        tangential displacement function = spot-weld-tangential
      end spot weld sw7
      #
      # Spot weld sw8
      begin spot weld sw8
        node set = nodelist_1018
        surface = surface_1007
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal
        tangential displacement function = spot-weld-tangential
      end spot weld sw8
      #
      # Spot weld sw9
      begin spot weld sw9
        node set = nodelist_1019
        surface = surface_1008
        failure decay cycles = 5
        failure envelope exponent = 2
        normal displacement function = spot-weld-normal
        tangential displacement function = spot-weld-tangential
      end spot weld sw9
      #
      # - NodeSet id 2 had name floor-bottom
      begin fixed displacement
        node set = nodelist_2
        components = X Y Z
      end fixed displacement
      #
      # - Block id 3 had name case
      begin initial velocity
        block = block_3
        direction = simba_v_shared_1
        magnitude = 10000
      end initial velocity
      #
      # - Block id 1 had name cover
      begin initial velocity
        block = block_1
        direction = simba_v_shared_1
        magnitude = 10000
      end initial velocity
      #
      # - Block id 11 had name 11
      begin initial velocity
        block = block_11
        direction = simba_v_shared_1
        magnitude = 10000
      end initial velocity
      #
      # - Block id 21 had name 21
      begin initial velocity
        block = block_21
        direction = simba_v_shared_1
        magnitude = 10000
      end initial velocity
      
      begin history output history-output
        database name = simple-assy-new-bolts.h
        database type = exodusII
        at step 0, increment = 20
        variable = node force_external at node 70299 as f_ext_sw1
        variable = node force_external at node 73267 as f_ext_sw2
        variable = node force_external at node 73627 as f_ext_sw3
        variable = node force_external at node 73447 as f_ext_sw4
        variable = node force_external at node 70659 as f_ext_sw5
        variable = node force_external at node 70539 as f_ext_sw6
        variable = node force_external at node 73867 as f_ext_sw7
        variable = node force_external at node 72967 as f_ext_sw8
        variable = node force_external at node 73087 as f_ext_sw9
        variable = node force_external at node 70119 as f_ext_sw10
        
      end history output history-output
      
    end presto region presto
    
  end presto procedure Apst_Procedure
  
end sierra simple-assy-new-bolts

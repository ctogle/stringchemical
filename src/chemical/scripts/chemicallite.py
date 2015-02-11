import modular_core.libfundamental as lfu
from modular_core.libfundamental import mobject as modular_object
from modular_core.libfundamental import run_parameter as run_param

import modular_core.libsimcomponents as lsc
import modular_core.libmath as lm
import modular_core.libgeometry as lgeo
import modular_core.libmodcomponents as lmc

import modular_core.fitting.libfitroutine as lfr
import modular_core.postprocessing.libpostprocess as lpp
import modular_core.criteria.libcriterion as lc

import stringchemical as chemfast
import stringchemical_timeout as chemfast_timeout

import sys, time, types, random, traceback
import numpy as np
from math import log as log
import cStringIO as sio

import pdb

if __name__ == 'chemical.scripts.chemicallite':
    lfu.check_gui_pack()
    lgm = lfu.gui_pack.lgm
    lgb = lfu.gui_pack.lgb
    lgd = lfu.gui_pack.lgd
if __name__ == '__main__':print 'stringchemical module'

module_name = 'chemical'
_system_string_ = '__initialized_system_string__'

class simulation_module(lmc.simulation_module):

    def _parse_variable(self,li,ensem,parser,procs,routs,targs):
        spl = lfu.msplit(li)
        name,value = spl
        varib = variable(name = name,value = value)
        return name,varib

    def _parse_function(self,li,ensem,parser,procs,routs,targs):
        spl = lfu.msplit(li)
        name,value = spl
        func = function(name = name,func_statement = value)
        return name,func

    def _parse_reaction(self,li,ensem,parser,procs,routs,targs):
        spl = lfu.msplit(li)
        rxn,lab = spl

        rxnspl = rxn.split(' ')
        arrows = ['<->','->','<-']
        for a in arrows:
            if a in rxnspl:
                divider = rxnspl.index(a)

        def stoch(side):
            sdxs = [num*2 for num in range(len(side)/2)]
            read = [(side[k + 1],int(side[k])) for k in sdxs]
            return read

        if rxnspl[divider] == '<->':
            r1,r2 = rxnspl[divider - 1],rxnspl[divider + 1]
            left  = [item for item in data[:divider-1] if not item == '+']
            right = [item for item in data[divider+2:] if not item == '+']
            left = stoch(left)
            right = stoch(right)
            rxn1 = reaction(name = lab+'1',rate = r1,used = left,produced = right)
            rxn2 = reaction(name = lab+'2',rate = r2,used = right,produced = left)
            return [rxn1,rxn2]

        elif rxnspl[divider] == '->':
            r1 = rxnspl[divider - 1]
            left  = [item for item in rxnspl[:divider-1] if not item == '+']
            right = [item for item in rxnspl[divider+1:] if not item == '+']
            left = stoch(left)
            right = stoch(right)
            rxn = reaction(name = lab,rate = r1,used = left,produced = right)
            return rxn

        elif data[divider] == '<-':
            r1 = rxnspl[divider + 1]
            left  = [item for item in rxnspl[:divider] if not item == '+']
            right = [item for item in rxnspl[divider+2:] if not item == '+']
            left = stoch(left)
            right = stoch(right)
            rxn = reaction(name = lab,rate = r1,used = right,produced = left)
            return rxn

    def _parse_species(self,li,ensem,parser,procs,routs,targs):
        spl = lfu.msplit(li)
        spec,value = spl
        new = species(name = spec,initial_count = value)
        return spec, new

    def __init__(self,*args,**kwargs):
        self.run_parameter_keys.extend(
            ['Variables','Functions','Reactions','Species'])
        self.parse_types.extend(
            ['variables','functions','reactions','species'])
        self.parse_funcs.extend(
            [self._parse_variable,self._parse_function, 
            self._parse_reaction,self._parse_species])
        lmc.simulation_module.__init__(self,*args,**kwargs)

    def _write_mcfg(self,mcfg_path,ensem):
        rparams = ensem.run_params
        mcfg = sio.StringIO()
        self._write_mcfg_run_param_key(rparams,'variables',mcfg)
        self._write_mcfg_run_param_key(rparams,'functions',mcfg)
        self._write_mcfg_run_param_key(rparams,'reactions',mcfg)
        self._write_mcfg_run_param_key(rparams,'species',mcfg)
        lmc.simulation_module._write_mcfg(mcfg_path,ensem,mcfg)



    # NEEDS SERIOUS REFACTORING
    # NEEDS SERIOUS REFACTORING
    # NEEDS SERIOUS REFACTORING
    # this is a handle to set param specific info per p-space location
    def _set_parameters(self):
        global _system_string_

        def make_rxn_string(rxn):
          used = '+'.join([''.join(['(', str(agent[1]), ')', 
              agent[0]]) for agent in rxn.used])
          prod = '+'.join([''.join(['(', str(agent[1]), ')', 
                            agent[0]]) for agent in rxn.produced])
          return '->'.join([used, str(rxn.rate), prod])

        def int_fix(cnt):
          if float(cnt) < 1: return 0
          else: return cnt

        params = ensem.run_params.partition['system']

        sub_spec = [':'.join([spec.label, str(int_fix(spec.initial_count))]) 
                          for spec in params['species'].values()]
        spec_string = '<species>' + ','.join(sub_spec)

        sub_var = [':'.join([key, str(var.value)]) for key, var in 
                                params['variables'].items()]
        variable_string = '<variables>' + ','.join(sub_var)

        def check_ext(afunc):
            extcnt = afunc.count('external_signal')
            fixed = []
            for exts in range(extcnt):
                leads = afunc.find('external_signal(')
                subfunc = afunc[leads+16:]
                presig = afunc[:leads+16]
                postsig = subfunc[subfunc.find('&'):]
                filename = subfunc[:subfunc.find(postsig)]

                #filename = subfunc[:subfunc.find('&')]
                with open(filename,'r') as handle:
                    extlines = handle.readlines()

                extstrx = sio.StringIO()
                extstry = sio.StringIO()
                for eline in extlines:
                    eline = eline.strip()
                    if not eline.count(',') > 0: continue
                    elx,ely = eline.split(',')
                    extstrx.write(str(elx))
                    extstrx.write('$')
                    extstry.write(str(ely))
                    extstry.write('$')

                fixhash = '%#%'
                extstrx.write('@')
                fixval = extstrx.getvalue() + extstry.getvalue()
                fixed.append((fixhash,fixval))
                afunc = presig + fixhash + postsig
            
            for fix in fixed: afunc = afunc.replace(fix[0],fix[1])
            #for fix in fixed: afunc = afunc.replace(fix[0],'---')
            return afunc
        
        sub_func = ['='.join([key, fu.func_statement.replace(',', '&')]) 
            for key, fu in params['functions'].items()]
        sub_func = [check_ext(sf) for sf in sub_func]
        function_string = '<functions>' + ','.join(sub_func)

        sub_rxn = ','.join([make_rxn_string(rxn) 
            for rxn in params['reactions']])
        reaction_string = '<reactions>' + sub_rxn

        sub_end = lc.read_criteria(params['end_criteria'], '')
        end_string = '<end>' + sub_end

        sub_capt = lc.read_criteria(params['capture_criteria'], '')
        capture_string = '<capture>' + sub_capt

        targs = params['plot_targets']
        sub_targ = ','.join(targs[3:] + targs[:3])
        target_string = '<targets>' + sub_targ + '||'

        system_string = spec_string + variable_string +\
            function_string + reaction_string + end_string +\
                capture_string + target_string

        _system_string_ = system_string
        



    def _reset_parameters(self):
        ensem = self.parent
        self._gui_memory()
        ensem.simulation_plan._reset_criteria_lists()
        ensem.run_params['variables'] = {}
        ensem.run_params['species'] = {}
        ensem.run_params['reactions'] = []
        ensem.run_params['functions'] = {}
        ensem.run_params['plot_targets'] = ['iteration','time']
        ensem.postprocess_plan._reset_process_list()
        output_plan = ensem.run_params['output_plans']['Simulation']
        output_plan.targeted = ['iteration','time']
        for w in output_plan.writers:w.targeted = ['iteration','time']

    def _gui_memory(self):
        self.module_memory = [
            lfu.data_container(selected_output_plan = 'Simulation', 
                selected_variable = 'None',selected_function = 'None', 
                selected_reaction = 'None',selected_species = 'None')]

    def _run_param_template(self,window,ensem,
            mobjname,key,handle_key,memory_key):
        new = (key,lgm.generate_add_remove_select_inspect_box_template(
            window = window,key = key,parent = ensem,
            labels = ['Add ' + mobjname,'Remove ' + mobjname], 
            wheres = [ensem.children,ensem.run_params[key]],
            selector_handle = (self.module_memory[0],handle_key),
            memory_handle = (self.module_memory[0],memory_key), 
            base_class = variable))
        return new

    def _panel_templates(self,*args,**kwargs):
        window = args[0]
        ensem = args[1]
        self._gui_memory()
        plot_target_labels = ['iteration','time'] +\
            ensem.run_params['species'].keys() +\
            ensem.run_params['variables'].keys() +\
            ensem.run_params['functions'].keys()
        panel_template_lookup =\
            lmc.simulation_module._panel_templates(
                self,window,ensem,plot_target_labels)
        panel_template_lookup.append(self._run_param_template(window,ensem,
            'Variable','variables','variable_selector','selected_variable'))
        panel_template_lookup.append(self._run_param_template(window,ensem,
            'Function','functions','function_selector','selected_function'))
        panel_template_lookup.append(self._run_param_template(window,ensem,
            'Reaction','reactions','reaction_selector','selected_reaction'))
        panel_template_lookup.append(self._run_param_template(window,ensem,
            'Species','species','species_selector','selected_species'))
        return panel_template_lookup

################################################################################
### run_parameter subclasses used to interface with the cython simulator
###   species
###   reaction
###   variable
###   function
################################################################################

class species(run_param):

    def __init__(self,*args,**kwargs):
        self._default('name','aspecies',**kwargs)
        self._default('initial_count',0,**kwargs)
        pspace_axes =\
          [lgeo.pspace_axis(instance = self,key = 'initial_count',
              bounds = [0,1000000],increment = 1,continuous = False)]
        self.pspace_axes = pspace_axes
        run_param.__init__(self,*args,**kwargs)

    def _string(self):
        return '\t' + self.name + ' : ' + str(self.initial_count)

    def _widget(self,*args,**kwargs):
        window = args[0]
        ensem = args[1]
        self._sanitize(*args,**kwargs)
        cartographer_support = lgm.cartographer_mason(window)
        self.widg_templates.append(
            lgm.interface_template_gui(
                mason = cartographer_support, 
                widgets = ['spin'], 
                instances = [[self]], 
                keys = [['initial_count']], 
                minimum_values = [[0]], 
                maximum_values = [[sys.maxint]], 
                initials = [[self.initial_count]], 
                box_labels = ['Initial Count'], 
                parameter_space_templates =\
                    [self.pspace_axes[0]]))
        self.widg_templates.append(
            lgm.interface_template_gui(
                keys = [['name']], 
                minimum_sizes = [[(150,50)]], 
                read_only = [True],
                instances = [[self]], 
                widgets = ['text'], 
                box_labels = ['Species Name']))
        run_param._widget(self,*args,from_sub = True)

################################################################################

class reaction(run_param):

    def __init__(self,*args,**kwargs):
        self._default('name','a reaction',**kwargs)
        self._default('rate',float(10.0),**kwargs)
        self._default('used',[],**kwargs)
        self._default('produced',[],**kwargs)
        pspace_axes =\
            [lgeo.pspace_axis(instance = self,key = 'rate',
                bounds = [0.0000000000001,100000000000.0], 
                continuous = True)]
        self.pspace_axes = pspace_axes
        run_param.__init__(self,*args,**kwargs)

    def _string(self):
        def _string_agents(agents):
            if not agents: return 'nothing'
            else: return ' '.join([str(a) for a in lfu.flatten(agents)])
        used_line = agents_to_line(self.used)
        produced_line = agents_to_line(self.produced)
        rxn_string = ' '.join([used_line,str(self.rate),'->',produced_line])
        rxn_string = '\t' + rxn_string + ' : ' + self.label
        return rxn_string

    def _widget(self,*args,**kwargs):
        window = args[0]
        ensem = args[1]
        spec_list = ensem.run_params['species'].keys()
        self.used = [u for u in self.used if u[0] in spec_list]
        self.produced = [p for p in self.produced if p[0] in spec_list]
        cartographer_support = lgm.cartographer_mason(window)
        self._sanitize(*args,**kwargs)
        left_template = lgm.interface_template_gui(
            panel_position = (0, 2), 
            mason = cartographer_support, 
            layout = 'vertical', 
            keys = [['name'],['rate']], 
            instances = [[self],[self]], 
            widgets = ['text','text'], 
            minimum_sizes = [[(400,100)],[(100,100)]], 
            box_labels = ['Reaction Name','Reaction Rate'], 
            initials = [[self.name],[self.rate]], 
            parameter_space_templates = [None,self.pspace_axes[0]])
        agents_template = lgm.interface_template_gui(
            panel_position = (0, 0), 
            layout = 'horizontal', 
            widgets = ['check_spin_list','check_spin_list'], 
            keys = [['used'],['produced']], 
            instances = [[self],[self]], 
            labels = [spec_list,spec_list],
            box_labels = ['Reagents','Products'])
        self.widg_templates.append(
            lgm.interface_template_gui(
                widgets = ['splitter'], 
                orientations = [['horizontal']], 
                templates = [[left_template,agents_template]]))
        run_param._widget(self,*args,from_sub = True)
        
################################################################################

class variable(run_param):
  
    def __init__(self,*args,**kwargs):
        self._default('name','a variable',**kwargs)
        self._default('value',1.0,**kwargs)
        pspace_axes = [
            lgeo.pspace_axis(instance = self,key = 'value',
                        bounds = [0.0,sys.float_info.max])]
        self.pspace_axes = pspace_axes
        run_param.__init__(self,*args,**kwargs)

    def _string(self):
        return '\t' + self.name + ' : ' + str(self.value)

    def _widget(self,*args,**kwargs):
        window = args[0]
        cartographer_support = lgm.cartographer_mason(window)
        self._sanitize(*args,**kwargs)
        self.widg_templates.append(
            lgm.interface_template_gui(
                widgets = ['spin'], 
                doubles = [[True]], 
                initials = [[float(self.value)]], 
                instances = [[self]], 
                keys = [['value']], 
                box_labels = ['Variable Value'], 
                mason = cartographer_support, 
                parameter_space_templates =\
                    [self.pspace_axes[0]]))
        self.widg_templates.append(
            lgm.interface_template_gui(
                widgets = ['text'], 
                read_only = [True],
                keys = [['name']], 
                instances = [[self]], 
                initials = [[self.name]], 
                box_labels = ['Variable Name']))
        run_param._widget(self,*args,from_sub = True)

################################################################################

class function(run_param):

    def __init__(self,*args,**kwargs):
        self._default('name','a function',**kwargs)
        self._default('func_statement','',**kwargs)
        run_param.__init__(self,*args,**kwargs)

    def _string(self):
        return '\t' + self.name + ' : ' + self.func_statement

    def _widget(self,*args,**kwargs):
        self._sanitize(*args,**kwargs)
        self.widg_templates.append(
            lgm.interface_template_gui(
                keys = [['func_statement']], 
                instances = [[self]], 
                widgets = ['text'], 
                minimum_sizes = [[(200,75)]], 
                box_labels = ['Function Statement'], 
                initials = [[self.func_statement]]))
        self.widg_templates.append(
            lgm.interface_template_gui(
                widgets = ['text'], 
                keys = [['name']], 
                instances = [[self]], 
                initials = [[self.name]], 
                box_labels = ['Function Name']))
        run_param._widget(self,*args,from_sub = True)

################################################################################
################################################################################





class sim_system(lsc.sim_system_external):

  def encode(self):
      global _system_string_
      self.system_string = _system_string_
      return

  def iterate(self):
    try:
      #seed = int(time.time()) + self.identifier
      if self.timeout:
        self.data = self.finalize_data(
          *chemfast_timeout.simulate(
            self.system_string, self.timeout))
            #self.system_string, seed, self.timeout))

      else:
        self.data = self.finalize_data(
          *chemfast.simulate(self.system_string))
          #*chemfast.simulate(self.system_string, seed))

    except ValueError:
      traceback.print_exc(file=sys.stdout)
      print 'simulation failed; aborting'
      self.bAbort = True
      raise ValueError
    except:
      traceback.print_exc(file=sys.stdout)
      print 'simulation failed; aborting'
      self.bAbort = True












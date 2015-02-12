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
import stringchemical_timeout as chemfast_to

import pdb,sys,time,types,random,traceback
import numpy as np
from math import log as log
from cStringIO import StringIO

if __name__ == 'chemical.scripts.chemicallite':
    lfu.check_gui_pack()
    lgm = lfu.gui_pack.lgm
    lgb = lfu.gui_pack.lgb
    lgd = lfu.gui_pack.lgd
if __name__ == '__main__':print 'stringchemical module'

module_name = 'chemical'

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
        new = species(name = spec,initial = value)
        return spec, new

    def __init__(self,*args,**kwargs):
        self._default('timeout',0.0,**kwargs)
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
        mcfg = StringIO()
        self._write_mcfg_run_param_key(rparams,'variables',mcfg)
        self._write_mcfg_run_param_key(rparams,'functions',mcfg)
        self._write_mcfg_run_param_key(rparams,'reactions',mcfg)
        self._write_mcfg_run_param_key(rparams,'species',mcfg)
        lmc.simulation_module._write_mcfg(mcfg_path,ensem,mcfg)

    # encode the species for cython simulation
    def _set_parameters_species(self,sysstr):
        species = self.parent.run_params['species']
        sub_spec = [species[s]._sim_string() for s in species]
        sysstr.write('<species>')
        sysstr.write(','.join(sub_spec))

    # encode the variables for cython simulation
    def _set_parameters_variables(self,sysstr):
        varis = self.parent.run_params['variables']
        sub_var = [varis[v]._sim_string() for v in varis]
        sysstr.write('<variables>')
        sysstr.write(','.join(sub_var))

    # encode the functions for cython simulation
    def _set_parameters_functions(self,sysstr):
        funcs = self.parent.run_params['functions']
        sub_fun = [funcs[f]._sim_string() for f in funcs]
        sysstr.write('<functions>')
        sysstr.write(','.join(sub_fun))

    # encode the reactions for cython simulation
    def _set_parameters_reactions(self,sysstr):
        rxns = self.parent.run_params['reactions']
        sub_rxs = [r._sim_string() for r in rxns]
        sysstr.write('<reactions>')
        sysstr.write(','.join(sub_rxs))

    # encode criteria for cython simulation
    def _set_parameters_criteria(self,sysstr,kind):
        crits = self.parent.run_params[kind+'_criteria']
        sub_crits = [c._sim_string() for c in crits]
        sysstr.write('<'+kind+'>')
        sysstr.write(','.join(sub_crits))

    # encode plot_targets for cython simulation
    def _set_parameters_plot_targets(self,sysstr):
        targs = self.parent.run_params['plot_targets']
        targs = targs[2:] + targs[:2]
        sysstr.write('<targets>')
        sysstr.write(','.join(targs))
        sysstr.write('||')

    # this is a handle to set param specific info per p-space location
    def _set_parameters(self):
        sysstr = StringIO()
        self._set_parameters_species(sysstr)
        self._set_parameters_variables(sysstr)
        self._set_parameters_functions(sysstr)
        self._set_parameters_reactions(sysstr)
        self._set_parameters_criteria(sysstr,'end')
        self._set_parameters_criteria(sysstr,'capture')
        self._set_parameters_plot_targets(sysstr)
        self.system_string = sysstr.getvalue()

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

    def _finalize_data(self,data,targets,ignore = []):
        reorder = []
        for name in self.parent.run_params['plot_targets']:
            if name in ignore or name not in targets: continue
            reorder.append(data[targets.index(name)])
        return np.array(reorder,dtype = np.float)

    def _simulate(self,*args,**kwargs):
        sstr = self.system_string
        if self.timeout:data = chemfast_to.simulate(sstr,self.timeout)
        else:data = chemfast.simulate(sstr)
        return self._finalize_data(*data)

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
        self._default('initial',0,**kwargs)
        pspace_axes =\
          [lgeo.pspace_axis(instance = self,key = 'initial',
              bounds = [0,1000000],increment = 1,continuous = False)]
        self.pspace_axes = pspace_axes
        run_param.__init__(self,*args,**kwargs)

    def _sim_string(self):
        return self.name + ':' + str(lfu.clamp(int(self.initial),0,sys.maxint))

    def _string(self):
        return '\t' + self.name + ' : ' + str(self.initial)

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
                keys = [['initial']], 
                minimum_values = [[0]], 
                maximum_values = [[sys.maxint]], 
                initials = [[self.initial]], 
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

    def _sim_string(self):
        def spec(agent):return '(' + str(agent[1]) + ')' + agent[0]
        def side(agents):return '+'.join([spec(a) for a in agents])
        return side(self.used)+'->'+str(self.rate)+'->'+side(self.produced)

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

    def _sim_string(self):
        return self.name + ':' + str(self.value)

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

    # modifies sim_string for ext_signal support
    #   THIS NEEDS MORE CLEANUP
    def _sim_string_ext_signal(self):
        afunc = self.func_statement
        extcnt = afunc.count('external_signal')
        fixed = []
        for exts in range(extcnt):
            leads = afunc.find('external_signal(')
            subfunc = afunc[leads+16:]
            presig = afunc[:leads+16]
            postsig = subfunc[subfunc.find('&'):]
            filename = subfunc[:subfunc.find(postsig)]
            with open(filename,'r') as handle:
                extlines = [l.strip() for l in handle.readlines()]

            extstrx = StringIO()
            extstry = StringIO()
            for eline in extlines:
                eline = eline.strip()
                if not eline.count(',') > 0: continue
                elx,ely = eline.split(',')
                extstrx.write(str(elx));extstrx.write('$')
                extstry.write(str(ely));extstry.write('$')

            fixhash = '%#%'
            extstrx.write('@')
            fixval = extstrx.getvalue() + extstry.getvalue()
            fixed.append((fixhash,fixval))
            afunc = presig + fixhash + postsig
            for fix in fixed: afunc = afunc.replace(fix[0],fix[1])
        return afunc

    def _sim_string(self):
        sysstr = self.name + '=' + self.func_statement.replaces(',','&')
        return self._sim_string_ext_signal(sysstr)

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











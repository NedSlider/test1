

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<ins> **Structure of the Program** </ins>

The program uses the *Agents* package
([*https://juliadynamics.github.io/Agents.jl/stable/*](https://juliadynamics.github.io/Agents.jl/stable/))
to simulate collections of cells. Agents provides a generic interface to
build ABM (Agent Based Modelling) applications. The main Agents routine,
**run!**(), calls two customised routines (**update_agent!**() and
**update_model!**()) to perform Agents Based Modelling simulations.

```julia
function run!(model,update_agent!,update_model!)
    n = 0          
    while ( n < maxCycles)
        for agent in allagents(model)            
        update_agent!(agent)            
    end                                                    
    update_model!(model)                      
    n += 1                                                
    end                                                                 
end
         Figure 1                                                                            
```

The basic Agents *run!*() function has been extended to allow users to
define an Agents model as a set of cells. Users control how cells behave
by attaching a set of 'Events' to each cell.


|aaa|bbb|ccc|

|a|b|c|



<table>
  <tr>
    <td align="centre" colspan="2">
    <pre lang="julia">
    function run!(model,update_cell!,update_model!)
        n = 0
        while ( n < maxCycles)
            for cell in allagents(model)
                update_cell!(agent)
            end 
        end
        update_model!(model)
        n += 1
    end
    </pre>
    </td>
    </tr>
   <tr>
   <td>
    <pre lang="julia">
    function update_cell!(agent, model)
        updateIntegration
        executeCellEvents(model,agent)
    end         
    </pre>
    </td>

   <td>
   <pre lang="julia">
    function() update_model!(model)
       model.nSteps += 1
      executeModelEvents(model)                         
    end
    </pre>                          
   </td>
  </tr>

  <tr>
  <td>
    <pre lang="julia">
    function executeCellEvents(model,cell)
        eventList = getfield(cell,:events_)
        for event in eventList
            if(event.test(model,cell,event))
                event.execute(model,cell,event)
            end
        end
    end
    </pre>
  </td>
  <td>
  <pre lang="julia">
   function executeModelEvents(model)
       cell = nothing
       eventList = model.events
       for event in eventList
           if(event.test(model,cell,event))
              event.execute(model,cell,event)
           end
       end
   end
  </pre>
  </td>
  </tr>
</table>


&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Figure 2

An 'Event' is a Julia struct consisting of a set of customised data and functions which
perform all of the required tasks for each Event.

```julia
mutable struct Event <: AbstractEvent
    name_::String                   # (unique) name of Event
    data_::Dict{String,Any}         # all data
    functions_::Dict{Symbol,Any}    # test(), execute(), save()
    global_::Bool                   # if true, 1 copy of the Event is shared between cells
end                                                     

Figure 2a     Event is a predefined Julia struct
```


Each **Event** has a unique name, a Dictionary (for storing all required
data) and a set of functions for performing tasks associated with the
Event. The 'global\_' variable is set to 'true' if an Event is shared
between several cells. If global\_ is false, each cell will receive it's
own copy of the Event. Global Events can be used to save memory when
performing generic tasks, e.g. writing variables to an output file.
There are three functions which must be defined:test(), execute(),
save(). By default, test() is assigned a dummy function (which always
returns true), and save() is assigned a dummy function which does
nothing. The execute() function <ins>must</ins> be defined. The test() function for each Event is always
executed. If test() returns true, the execute() function is called. The
save() function is used to save data (e.g. write variables to a file).
Events provide a simple mechanism which allows users to customise a
simulation by choosing different events, Figure 2. This mechanism also
allows programmers to add functionality by creating new Events.

<table>
  <tr>
  <td>
    <pre lang="julia">
    function timToDivide(model,cell,event)
    eventData = getfield(event,:data_)
       nutrient = Events.getEventVariable("nutrient",eventData)
       nCells = length(model.agents)
       max = 0.05 - nCells * nutrient
       probability = nSteps * rand(Uniform(0.0,max))
       if(probability > 1.0)                        
           return(true)
       end
       return(false)
    end
    </pre>
  </td>
  <td>
      <pre lang="julia">
      function divide_cell(model,cell,event)
          model.lastCell += 2
          data = getfield(event,:data_)
          fraction = Events.getEventVariable("fraction",data)
          u1 = divideResources(cell,fraction)
          u2 = divideResources(cell,1.0-fraction)
          cellLineage = getfield(cell,:lineage_)
          cellLineage.status_ = "divided"
          cell1 = createNewCell(model,cell,nextAgent,u1,(0.0,0.0),cellIndex)
          cell2 = createNewCell(model,cell,nextAgent+1,u2,(0.0,0.0),cellIndex)
          kill_agent!(cell,model)
      end
      </pre>
  </td>
  </tr>
</table>

An example of two functions, timeToDivide() and divide_cell() which can be used to create an Event which decides when a cell should divide and then creates two new daughter cells. The variables, nutrient and fraction, are stored in the Event data Dict. They are updated, when necessary, by the Event.

  **Main Input File**


```
  Cell:A cell_a.dat                                            # Catalyst input data
  CellEvents:cell_events.dat                                   # contains all predefined cell events
  ModelEvents:model_events.dat                                 # contains all predefined model events
  Events:A updateGrowthA{updateGrowth(growthLimit=0.01,fraction=0.55)}
  Events:model updateNutrientA{updateNutrient()}
  Events:model saveCellData{saveData()}
  Events:model cellCount{modelCellCount()}                     # events attached to model
  Initialise:demo.cells.dat                                    # contains the initial cells
  Output:./demo file-prefix
  u0:A steadystate.csv 1.0
  Write:Cell A r e_t e_m q m_r m_t m_m m_q c_r c_t c_m c_q a s_i
```

The main input file uses a set of keywords:**Cell**, **CellEvents**, **ModelEvents**, **Events**, **Initialise**,
**u0** , **Write**, to define all of the required data.

**1) Cell**
Defines ODEs. There are two mechanisms for creating ODEs

## ***1a Define ODEs using Catalyst Reaction Networks***

**Cell**:A cell_cat.dat

This command defines the input for creating ODEs using Catalyst

Format: Cell:\<*cell-label*\> \<*file*\>

Cell should be followed immediately by a colon and a Cell 'label', and then a data file which defines the reaction network in Catalyst format (see <https://github.com/SciML/Catalyst.jl>)

The format for the Catalyst data is defined in [Appendix 1a](#appendix-1a) .

The cell *label* must be unique because it allows users to define multiple cell types which have different Catalyst reactions.

## ****1b Defining Customised ODE functions****

**Cell**:A cell_ode.dat

This command defines the input for creating ODEs using customised Julia functions.

The functions must be defined in the file odeFunctions.jl

The format for ODEs functions is define in [Appendix 1b](#appendix-1b).

## 2 CellEvents

**CellEvents**:cell_events.dat

This command allows a user to define a set of Events which can be
attached to any cell.

Format *CellEvents*:\<*filename*\>

The key word, *CellEvents*, should be followed immediately by a colon,
and the name of the file containing the definitions for each CellEvent.
The format of this file is defined in [Appendix 2](#appendix-2)

The data and functions defined in the 'cell events' file are used to create cell events dynamically at run time. The Events are then stored in a Dict.

## **3 ModelEvents**

**ModelEvents**:model_events.dat

This command allows users to define a list of all required model Events.

Format *ModelEvents:*\<f*ilename*\>

The key word, *ModelEvents*, must be followed by a colon, and then the name of the file which contains definitions for a set of ModelEvents. The format of this file is defined in [Appendix 3](#appendix-3) .

The data and functions defined in the 'model events' file are used to create model events dynamically at run time. The Events are then stored in a Dict.

**4. Creating Instances of Cell Events**

An Event *Instance* is an Event which has been created by copying a
predefined Event from the Events Dict and then modifying it by resetting
some of the variables stored in the Event data Dict. This allows users
to define customised Events at run time.

The 'Events' keyword is used to create instances of Events.

**Events**:A updateGrowthA{updateGrowth(growthLimit=1.0,fraction=0.65)}

Format Events:\<*cell label*\>
\<*instance-name*\>**{**\<*event-name*\>**(**comma separated
variables**)}**

The Events keyword allows users to assign specific instances of Events
to a set of cells with the same \<cell label\>.

This keyword is followed by a colon and then a cell label (e.g. A).

In the example above the \<*instance-name*\> is given a (unique)
name:*updateGrowthA*.

The remainder of the input enclosed by brackets, **{}**, defines an
instance of the Event, e.g.

{updateGrowth(growthLimit=1.0,fraction=0.65)}

This must correspond to a predefined Event (e.g. updateGrowth,
previously defined in the CellEvents file; see above).

This is subsequently parsed by removing the brackets, {}.

The remainder,

updateGrowth(growthLimit=1.0,fraction=0.65)

defines a Cell Event name, *updateGrowth*,

and a set of customised values for the Event variables:
(*growthLimit*=1.0,*fraction*=0.65).

The variable names and types must also be consistent with the original
definition of the Event.

An instance of this Event will be created at run time and assigned to
those cells which have the correct \<*cell-label*\>

**5. Creating Instances of Model Events**

Instances of Model Events can also be created using the Events keyword

Events:model update_nutrient{updateNutrient()}.

The only difference between this and the previous 'Cell Event' example
is the cell label has been replaced by 'model'. In this example, the
program will create an instance (labelled update_nutrient) of the
updateNutrient Event at runtime and attach it to the model.

**6. Initialising Cells**

Initialise:demo2.cells.dat

This command is used to define the initial set of cells.

Format *Initialise*:\<*cell-file*\>

\<*cell-file*\> contains the initial set of cells.

The format of the cell file is:

Cell \<*cell-label*\> \<*cell-id*\> (xccord, ycoord).

\<*cell-id*\> is a unique integer

e.g. Cell A 1 (x,y)

Cells defined in this file are used to build the ABM model. First a
model is created. Next the cells defined in this file are added. Each
cell is given a copy of the relevant CellEvent instances. Then the
ModelEvent instances are attached to the model. A serialized model is
then written to a .jls file.

**7. Defining Output Variables**

The **Write** keyword is used to define a list of variables to output,
e.g.

Write:Cell A rmr em rmq rmt et rmm mt mm q si mq mr r a rep mrep rmrep

Format Write:Cell \<*cell-label*\>
\<*list-of-space-separated-variables*\>

## **Appendix 1a**

[return to main text](#1a-Define-odes-using-catalyst-reaction-networks)



Example input file for Catalyst reactions

```
catalyst:catalyst.dat
p       s=1e4  d_m=0.1 n_s=0.5 n_r=7459  n_t=300 n_m=300 n_q=300 γ_max=1260  K_γ=7   v_t=726 K_t=1000    v_m=5800 K_m=1000        w_r_max=930     w_t_max=4.14    w_m_max=4.14    w_q_max=949 θ_r=427 θ_t=4.38        θ_m=4.38        θ_q=4.38        K_q=152219      h_q=4  kb=0.95e-2      ku=1    M=1e8
u0      r=10   e_t=0   e_m=0   q=0    m_r=0   m_t=0   m_m=0  m_q=0   c_r=0   c_t=0  c_m=0   c_q=0   a=1000 s_i=0
```

The input file uses three keywords (**catalyst:,p,u0**)

**p** defines a set of p values

p is immediately followed by a space and then a set of p-values:
\<*name*\>=\<*value*\> (no spaces)

The defined p-values are separated by spaces.

**u0** defines initial values for u0, the values are defined as above.

**catalyst:** specifies the set of reactions (as specified in the
Catalyst documentation)

Reaction Network defined in:catalyst.dat

```
rn_lamda = @reaction_network begin
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), r --> ∅
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), e_t --> ∅
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), e_m --> ∅
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), q --> ∅
        ((w_r_max * a)/(θ_r + a)), ∅ --> m_r
        ((w_t_max * a)/(θ_t + a)), ∅ --> m_t
        ((w_m_max * a)/(θ_m + a)), ∅ --> m_m
        ((w_q_max * a)/((θ_q + a)*(1 + (q/K_q)^h_q))), ∅ --> m_q
        (((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M) + d_m), m_r --> ∅
        (((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M) + d_m), m_t --> ∅
        (((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M) + d_m), m_m --> ∅
        (((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M) + d_m), m_q --> ∅
        (kb,ku), m_r + r <--> c_r
        (kb,ku), m_t + r <--> c_t
        (kb,ku), m_m + r <--> c_m
        (kb,ku), m_q + r <--> c_q
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), c_r --> ∅
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), c_t --> ∅
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), c_m --> ∅
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), c_q --> ∅
        (γ_max/(n_r * (K_γ + a))), a + c_r --> r + m_r + r
        (γ_max/(n_t * (K_γ + a))), a + c_t --> r + m_t + e_t
        (γ_max/(n_m * (K_γ + a))), a + c_m --> r + m_m + e_m
        (γ_max/(n_q * (K_γ + a))), a + c_q --> r + m_q + q
        ((γ_max * c_r * (n_r - 1))/(n_r * (K_γ + a))), a --> 0
        ((γ_max * c_t * (n_t - 1))/(n_t * (K_γ + a))), a --> 0
        ((γ_max * c_m * (n_m - 1))/(n_m * (K_γ + a))), a --> 0
        ((γ_max * c_q * (n_q - 1))/(n_q * (K_γ + a))), a --> 0
        ((e_m*v_m)/(K_m + s_i)), s_i --> n_s*a
        ((e_t*v_t*s)/(K_t + s)), 0 --> s_i
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), s_i --> ∅
        ((a*γ_max/((K_γ + a)))*(c_r + c_t + c_m + c_q)/M), a --> ∅
end
```

## Appendix 1b
[return to main text](#1b-defining-customised-ode-functions)

Example input for defining customised Julia ODE functions.


The keyword, *odeFunction*, (followed by colon) specifies a function
defining ODEs (e.g. *testFunction!*) which must be defined in
odeFunctions.jl

```
odeFunction:testFunction!
p	thetar=426.8693338968694	k_cm=0.005990373118888	s0=10000.0	
gmax=1260.0	cl=0.0	thetax=4.379733394834643	Kt=1000.0	M=1.0e8
we=4.139172187824451	Km=1000.0	vm=5800.0	nx=300.0	Kq=15221
9.040373749	Kp=180.1378030928276	vt=726.0	wr=929.9678874564831	
wq=948.9349882947897	wp=0.0	nq=4.0	nr=7549.0	ns=0.5	kurep=10000.0	
kbrep=1.0	wmaxrep=1.0	dmrep=0.34657359027997264	dprep=0.17328679
513998632	
u0	rmr=807.530517658162	em=7066.62814403594	rmq=2316.16746788774	
rmt=69.1538861165317	et=7066.62758622631	rmm=69.1538892849256	mt=9.640
96853393699	mm=9.64097064704617	q=236681.293193619	si=128.404551112
062	mq=322.904581569518	mr=5.94267303259607	r=17.3282796528522	
a=9.17329183561663	rep=0.0	mrep=0.0	rmrep=0.0
```
## **Appendix 2**     

[return to main text](#2-cellevents)

Example ‘Cell’ Events file

```
Events

CellEvent:updateGrowth
data:growthLimit(Float64) = 1.0
data:fraction(Float64) = 0.60
cell:growthProg(Float64) = 0.0
data:growthIncrease(Float64) = 0.0
reset:growthProg = 0.0, fraction = 0.5
equation:growthRateEqn = ((c_q + c_m + c_t + c_r) * (γ_max*a/(K_γ + a)))/M
test:checkGrowth
execute:divideByGrowth
save:dontSave

end
```

The first line is always *Events*, and the last line is always *end*.

The input for each is Event is controlled by several keywords: **CellEvent**, **data**, **cell**, **model**, **equation**, **test**, **execute**, **save**.

Keywords are always followed by a colon, e.g. **CellEvent**:updateGrowth

Format *CellEvent*:\<*event-name*\>

*CellEvent* is followed by a colon and then a name for the Event. The
name must be unique.

**data**:growthLimit(Float64) = 1.0

The data keyword is used to define 'Event' data. Event data are stored in the Event data dictionary.

Format data:\<*variable name*\>**(**\<*variable type*\>**)** =
\<*value*\>

\<*variable name*\> must be unique and \<*value*\> must correspond to a
valid \<*variable type*\>

**cell**:growthProg(Float64) = 0.0

Format cell:\<*variable name*\>**(**\<*variable type*\>**)** =
\<*value*\>

The cell keyword is similar to the data keyword; the format is
identical.

The input is parsed and a new variable is attached to the cell.

This allows users to create their own cell variables.

**equation**:growthRateEqn = ((c_q + c_m + c_t + c_r) \*
(γ_max\*a/(K_γ + a)))/M

The *equation* keyword allows users to create their own equations, which
are then accessible to any function associated with the Event.

Format *equation*:\<*equation name*\> **=**\<*equation definition*\>

Equations are parsed and stored in the Event data dictionary.

The equation name must correspond to an equation required by at least
one of the Event functions. Functions which use this equation will
retrieve it, by name, from the data dictionary and evaluate it.

**reset**:growthProg=0.0, fraction=0.5

The *reset*: keyword is used to re-initialise events. This is most
likely to happen during cell division. When a cell divides each daughter
cell gets a copy of the parent's events. Some of these events may
require reinitialisation. In this example a CellEvent is monitoring the
growth of the current (parent) cell and the value of *growthProg* is
constantly changing. When the parent cell divides each daughter cell
inherits a list of Events from the parent cell. Many of the Event
variables should be reset. The default reset values for all the Event
variables are the initial values given to each variable. The reset:
keyword allows users to specify alternative values. In this example the
reset: keyword is used to define new initial values for *growthProg* and
*fraction*.

reset:growthProg = 0.0,fraction=0.5

Now, when the parent cell divides, *growthProg* will be set to zero and
*fraction* to 0.7

Format *reset*: \<v*ariable name 1*\> = \<*value 1*\> , \<*variable name
2*\> = \<*value 2*\>

Variables are always separated by commas.

**test**:checkGrowth

Format *test*:\<f*unction name*\>

test, defines the function an Event will use as the test function. The function name must correspond to an existing function. By default, a dummy function which always returns 'true' is used if test has not been
defined.

**execute**:divideByGrowth

Format *execute*:\<*function name*\>

execute, defines the function an Event will execute (if test() is true).

**save**:dontSave

Format save:\<*function name*\>

save, defines the function an Event will use to save data. By default, a dummy function which does nothing is used if save has not been defined.

## **Appendix 3**

[return to main text](#3-modelevents)

Example ‘Model’ Events file


```
Events

 ModelEvent:updateNutrient
 data:cellCount(Int32) = 1
 model:nutrient(Float64) = 1000.0
 data:totalProduct(Float64) = 0.0
 data:productSymbol(Symbol) = s_i
 execute:updateNutrient

 end
```

The first line is always 'Events' and the last line is always 'end'. The
format of this file is the same as appendix 2. The same keywords are
used. This example also uses the 'model' keyword.

**model**:nutrient(Float64) = 1000.0

This creates a new variable, *nutrient*, which is attached to the model.
In the above example the execute function is assigned to the function,
*updateNutrient*().




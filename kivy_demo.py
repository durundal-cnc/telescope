from kivy.app import App, async_runTouchApp
from kivy.uix.gridlayout import GridLayout
from kivy.uix.label import Label
from kivy.uix.textinput import TextInput
from kivy.uix.button import Button
from kivy.uix.image import AsyncImage
from kivy.uix.dropdown import DropDown
from kivy.uix.checkbox import CheckBox
from kivy.metrics import dp
from kivy.uix.boxlayout import BoxLayout
from kivy.uix.progressbar import ProgressBar
from kivy.uix.checkbox import CheckBox

from kivy.clock import Clock


from functools import partial
import trio #for asynch running of the telescope thread



from kivy.properties import (
    NumericProperty,
    StringProperty,
    BooleanProperty,
    ListProperty,
)
import re
import sys


sys.path.append(r'/Users/andrewmiller/telescope/')

#my functions
import config # https://docs.python.org/3/faq/programming.html#how-do-i-share-global-variables-across-modules



config.x = 'myvar'

async def run_app_GUI(root, nursery):
    '''This method, which runs Kivy, is run by trio as one of the coroutines.
    '''
    # trio needs to be set so that it'll be used for the event loop
    await async_runTouchApp(root, async_lib='trio')  # run Kivy, asynch_runTouchApp is kivy function
    print('App done')
    # now cancel all the other tasks that may be running
    nursery.cancel_scope.cancel()



async def run_telescope(root):
    '''This method is also run by trio and periodically prints something.'''
    try:
        while True:
            print('Sitting on the beach')
            root.console.text += 'Sitting on the beach'

            config.x = 'reset'
            print(config.x)
            root.console.text += config.x

            await trio.sleep(2)
            x = 0
            while x<100000:
                x = x+1
                # print(x)
                # print(config.x)
                root.iPhone_stats.text = str(x)
                root.roboclaw_stats.text = str(x)
                
                root.trackprogressbar.value = x
                root.trackprogressbar.max = 100000
                await trio.sleep(0) #checkpoint without blocking (e.g. GUI can operate now)
    except trio.Cancelled as e:
        print('Wasting time was canceled', e)
    finally:
        # when canceled, print that it finished
        print('Done wasting time')

# async def read_log(root, nursery):
#     '''This method, which runs Kivy, is run by trio as one of the coroutines.
#     '''
#     try:
#         while True:
#             print('Reading console log')
#             await trio.sleep(2)
            
#             root.console.text = 

#             await trio.sleep(0) #checkpoint without blocking (e.g. GUI can operate now)
#     except trio.Cancelled as e:
#         print('Console log cancelled', e)
#     finally:
#         # when canceled, print that it finished
#         print('Done with console log')



class FloatInput(TextInput): #text input that only accepts numbers

    pat = re.compile('[^0-9]')
    def insert_text(self, substring, from_undo=False):
        pat = self.pat
        if '.' in self.text:
            s = re.sub(pat, '', substring)
        else:
            s = '.'.join(
                re.sub(pat, '', s)
                for s in substring.split('.', 1)
            )
        return super().insert_text(s, from_undo=from_undo)

class MainScreen(GridLayout):

    


    def __init__(self, **kwargs):
        super(MainScreen, self).__init__(**kwargs)
        self.cols = 3
        self.add_widget(Label(text='User Name'))
        self.username = TextInput(multiline=False)
        self.add_widget(self.username)
        self.add_widget(Label(text='password'))
        self.password = TextInput(password=True, multiline=False)
        self.add_widget(self.password)
        
        target_pairs = []

        
        #input_type is an OptionsProperty and defaults to ‘null’. Can be one of ‘null’, ‘text’, ‘number’, ‘url’, ‘mail’, ‘datetime’, ‘tel’ or ‘address’.
        #TextInput(hint_text='Kgs', input_filter = 'float', multiline=False, write_tab=False) #numeric only, tab moves to next object instead of writing \tab

# def on_focus(instance, value):
#     if value:
#         print('User focused', instance)
#     else:
#         print('User defocused', instance)

# textinput = TextInput()
# textinput.bind(focus=on_focus)


# def on_text(instance, value):
#     print('The widget', instance, 'have:', value)

# textinput = TextInput()
# textinput.bind(text=on_text)

        def manual_Az_focus_change(instance, value):
            print(config.x)
            self.manual_Az.text = str( float(self.manual_Az.text) %360)
            config.x = 'manual az config.x'

        self.add_widget(Label(text='Manual Az (deg)'))
        self.manual_Az = TextInput(multiline=False, hint_text='Manual Az (deg)', input_filter = 'float', write_tab = False) #input for numeric manual Az degrees
        self.manual_Az.text = '0.0'
        #self.manual_Az.bind(text = manual_Az_inputchange) #executes callback any time text is changed
        self.manual_Az.bind(focus = manual_Az_focus_change) #executes callback any time focus is changed
        self.add_widget(self.manual_Az)


        def manual_El_focus_change(instance, value):
            print(config.x)
            self.manual_El.text = str( float(self.manual_El.text) %360)
            config.x = 'manual el config.x'
        self.add_widget(Label(text='Manual El (deg)'))
        self.manual_El = TextInput(multiline=False, hint_text='Manual El (deg)', input_filter = 'float', write_tab = False) #input for numeric manual El degrees
        self.manual_El.text = '0.0'
        self.manual_El.bind(focus = manual_El_focus_change)
        self.add_widget(self.manual_El)

        def on_manual_AzEl_enter(instance): #button to accept manual typed Az/El
            #print('User pressed enter in', instance)
            #coerce output to 0-360
            self.manual_Az.text = str( float(self.manual_Az.text) %360)
            self.manual_El.text = str( float(self.manual_El.text) %360)

            print('Az : ' + self.manual_Az.text)
            print('El : ' + self.manual_El.text)
        self.manual_AzEl = (Button(text='Move AzEl'))
        self.manual_AzEl.bind(on_press=on_manual_AzEl_enter)
        self.add_widget(self.manual_AzEl)

#self.mytext.bind(text = self.calc)


    #     #Clock.schedule_interval(self.update, 1) #can use for timing updates
    #     return self.layout

    # def update(self, *args):
    #     self.name.text = str(self.current_i)
    #     self.current_i += 1
    #     if self.current_i >= 50:
    #         Clock.unschedule(self.update)

        #manual control buttons
        def plus_Az_callback(instance):
            self.manual_Az.text = str(float(self.manual_Az.text) + 1)
            print('The button <%s> is being pressed' % instance.text)
            self.manual_AzEl.trigger_action(0.1) # argument is for how long button should be pressed

        def minus_Az_callback(instance):
            self.manual_Az.text = str(float(self.manual_Az.text) - 1)
            print('The button <%s> is being pressed' % instance.text)
            self.manual_AzEl.trigger_action(0.1) # argument is for how long button should be pressed

        def plus_El_callback(instance):
            self.manual_El.text = str(float(self.manual_El.text) + 1)
            print('The button <%s> is being pressed' % instance.text)
            self.manual_AzEl.trigger_action(0.1) # argument is for how long button should be pressed

        def minus_El_callback(instance):
            self.manual_El.text = str(float(self.manual_El.text) - 1)
            print('The button <%s> is being pressed' % instance.text)
            self.manual_AzEl.trigger_action(0.1) # argument is for how long button should be pressed

        self.plus_Az_button = Button(text='+Az')
        self.plus_Az_button.bind(on_press=plus_Az_callback)
        self.add_widget(self.plus_Az_button)

        self.minus_Az_button = Button(text='-Az')
        self.minus_Az_button.bind(on_press=minus_Az_callback)
        self.add_widget(self.minus_Az_button)

        self.plus_El_button = Button(text='+El')
        self.plus_El_button.bind(on_press=plus_El_callback)
        self.add_widget(self.plus_El_button)

        self.minus_El_button = Button(text='-El')
        self.minus_El_button.bind(on_press=minus_El_callback)
        self.add_widget(self.minus_El_button)

        ImageUrl = 'https://kivy.org/doc/stable/_static/logo-kivy.png'
        self.add_widget(AsyncImage(source=ImageUrl)) #display camera images

        self.iPhone_stats = Label(text='iPhone boot text')
        self.add_widget(self.iPhone_stats)
        self.roboclaw_stats = Label(text='roboclaw boot text')
        self.add_widget(self.roboclaw_stats)
        self.az_PV = Label(text='az PV boot text')
        self.add_widget(self.az_PV)
        self.el_PV = Label(text='el PV boot text')
        self.add_widget(self.el_PV)
        self.current_target_stats = Label(text='Current target:')
        self.add_widget(self.current_target_stats)

        # def target_list(instance, value):
        #     print(config.x)
        #     self.target_list.text = str( float(self.manual_El.text) %360)
        #     config.x = 'manual el config.x'
        self.add_widget(Label(text='Target list (name, method)'))
        self.target_list_textinput = TextInput(multiline=False, hint_text='List of targets', write_tab = False) #input for numeric manual El degrees
        self.target_list_textinput.text = 'sun astronomy moon astronomy'
        self.add_widget(self.target_list_textinput)
        
        self.console = Label(text='Console output')
        self.add_widget(self.console)
        
###works up until here (dropdown creation)

        self.dropdown = DropDown()
        def dropdown_callback(instance, value):
            print('inside dropdown callback')
            setattr(self.choose_target, 'text', value)
            print('value' + value)
        self.dropdown.bind(on_select=dropdown_callback)



        def compute_targets_callback(instance):
            print('computing target coorindates')
            target_list_str = self.target_list_textinput.text
            print(target_list_str)
            target_list = target_list_str.split() #= 'moon astronomy polaris astronomy sun astronomy ISS satellite_tracking'
            print(target_list)
            target_pairs = []
            a = iter(target_list)
            for x,y in zip(a, a):
                target_pairs.append([x,y])
            print(target_pairs)
#            time_of_individual_tracks, name_of_individual_tracks, individual_tracks = compute_target_coords(target_list = target_pairs)
            #add to dropdown target selector
            
            ######pass target_pairs to target computing function here
            
            #create a dropdown with 10 buttons
            print('len of target pairs')
            print(len(target_pairs))
            self.dropdown.clear_widgets() #remove the old target widgets from the dropdown

            #create a dropdown with 10 buttons
            for index in range(len(target_pairs)):
                # When adding widgets, we need to specify the height manually
                # (disabling the size_hint_y) so the dropdown can calculate
                # the area it needs.
            
                target = Button(text='Value %d' % index + target_pairs[index][0] , size_hint_y=None, height=44)
            
                # for each button, attach a callback that will call the select() method
                # on the dropdown. We'll pass the text of the button as the data of the
                # selection.
                target.bind(on_release=lambda btn: self.dropdown.select(btn.text))
            
                # then add the button inside the dropdown
                self.dropdown.add_widget(target)
            
            
        self.compute_targets = Button(text='Compute target coords', size_hint_y = None)
        self.compute_targets.bind(on_press=compute_targets_callback)
        self.add_widget(self.compute_targets)


        def halt_telescope_callback(instance):
            print('halting telescope')
            #AzSV = AzPV
            #ElSV = ElPV
            #mode = 'manual'
        self.halt_telescope = Button(text='Halt telescope', size_hint_y = None)
        self.halt_telescope.bind(on_press=halt_telescope_callback)
        self.add_widget(self.halt_telescope)

        #choose new target


        def choose_target_callback(instance):
            
            #create a dropdown with 10 buttons
            print('len of target pairs')
            print(len(target_pairs))
                

        self.choose_target = Button(text='Select target', size_hint_y= None) #this is the main button for the dropdown (which contains other buttons and is hidden)
        self.choose_target.bind(on_release=self.dropdown.open)
        self.add_widget(self.choose_target)

        self.trackprogressbar = ProgressBar(max=1000)#max=len(THE LENGTH OF THE TRACK))
        
        self.trackprogressbar.value = 750
        self.trackprogressbar.max = 10000
        self.add_widget(self.trackprogressbar)
        
    

        def on_auto_checkbox_active(checkbox, value):
            if value:
                print('The checkbox', self.automanual_checkbox, 'is active')
            else:
                print('The checkbox', self.automanual_checkbox, 'is inactive')

        self.automanual_checkbox = CheckBox()
        self.automanual_checkbox.bind(active=on_auto_checkbox_active)
        self.add_widget(self.automanual_checkbox)

        self.automanual_label = Label(text='Auto next target')
        self.add_widget(self.automanual_label)

# class MyApp(App):
#     print('Doing some stuff')
#     def build(self):
#         tele = MainScreen()    
#         tele.manual_El.text = '123'

#         return tele
    
# if __name__ == '__main__':
#     MyApp().run()
    
        

if __name__ == '__main__':
    async def root_func():
        '''trio needs to run a function, so this is it. '''

        root = MainScreen()    #Builder.load_string(kv)  # root widget
        async with trio.open_nursery() as nursery:
            '''In trio you create a nursery, in which you schedule async
            functions to be run by the nursery simultaneously as tasks.

            This will run all two methods starting in random order
            asynchronously and then block until they are finished or canceled
            at the `with` level. '''
            nursery.start_soon(run_app_GUI, root, nursery)
            nursery.start_soon(run_telescope, root)

    trio.run(root_func)
    
    
    
    
    



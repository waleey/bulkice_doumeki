import numpy as np
from Event import Event

class NUANCEReader:
    """Simple Reader for NUANCE format
    """
    def __init__(self, file_path):   #this method will always execute everytime an instance is created
        self.file_path = file_path
        self.file_object = None
        self.n_read = 0
        self.events = []
        pass

    def __enter__(self):   #this method will always execute everytime an instance is created
        self.file_object = open(self.file_path, 'r')
        return self

    def __exit__(self, type, value, traceback):  #this method will always execute everytime an instance is created
        # Add exception handling here, as needed
        self.file_object.close()

    def get_events(self):
        """Get Events from Open NUANCE file
        """
        events = []
        lines = []

        for line in self.file_object.readlines():
            if line[0] == '#':
                continue
            elif 'begin' in line:
                evt_lines = []  # start new event lines
                continue
            elif 'end' in line:
                events.append(Event.from_text(evt_lines))
                continue
            elif 'stop' in line:
                continue
            evt_lines.append(line.strip()[2:])  # Removes leading '$ ' and trailing '\n'
        return events

import os
import subprocess
from abc import ABCMeta, abstractmethod
from pathlib import Path


class Tool:

    __metaclass__ = ABCMeta

    def __init__(self, exe: str, full_name: str, exe_path: str):
        self.exe = exe
        self.full_name = full_name
        self.exe_path = exe_path

    def run(self, parameters: dict):
        if self.__exe_exists():
            command = self.__create_command(parameters)
            return self.run_command(command)
        else:
            raise Exception('{0} executable could been found on path {1}'.format(self.exe, self.exe_path))

    @abstractmethod
    def run_command(self, command):
        pass

    def __create_command(self, parameters: dict) -> str:
        return self.exe + ''.join(' {} {} '.format(key, val) for key, val in parameters.items())

    def __exe_exists(self) -> bool:
        return Path(self.exe_path).is_file()


class StrikeEx(Tool):

    def __init__(self, exe: str = 'strike', full_name: str = 'Single structure induced evaluation',
                 exe_path: str = '/usr/local/bin/strike'):
        super(StrikeEx, self).__init__(exe, full_name, exe_path)

    def run_command(self, command) -> float:
        bytess = subprocess.check_output(command, shell=True, env=os.environ.copy())
        return float("".join(map(chr, bytess)).split('\n')[-2])

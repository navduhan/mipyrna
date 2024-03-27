#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
:Title: This is loggin module for miPyRNA

:Author : Naveen Duhan
'''


import logging

def MiPyRNALogger(mode, log):

    """ This function intialize logger in the miPyRNA modules

    :param  mode: Logger name for the module

    :param  log: File name for logging
    """
    logger = logging.getLogger(log)
    logger.propagate=False

    # set format for logging
    logFormatter =logging.Formatter('[%(asctime)s]  %(module)s :: %(levelname)s : %(message)s',datefmt='%H:%M:%S')
    
    logger.setLevel(logging.DEBUG)
    # write log in a user provided log file
    
    fileHandler = logging.FileHandler("{}".format('mipyrna'), mode= mode)

    fileHandler.setFormatter(logFormatter)

    logger.addHandler(fileHandler)

    consoleHandler = logging.StreamHandler()

    consoleHandler.setFormatter(logging.Formatter('[%(asctime)s]  %(module)s :: %(levelname)s : %(message)s',datefmt='%H:%M:%S'))

    consoleHandler.setLevel(logging.DEBUG)

    logger.addHandler(consoleHandler)

    return logger
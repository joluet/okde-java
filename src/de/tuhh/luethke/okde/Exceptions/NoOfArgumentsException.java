/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.Exceptions;

public class NoOfArgumentsException extends Exception {

	private static final long serialVersionUID = -790450771519569118L;

	public NoOfArgumentsException() {
	}

	public NoOfArgumentsException(String msg) {
		super(msg);
	}
}

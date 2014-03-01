/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.Exceptions;

public class TooManyComponentsException extends Exception {

	private static final long serialVersionUID = 5182689601412391328L;

	public TooManyComponentsException() {
	}

	public TooManyComponentsException(String msg) {
		super(msg);
	}
}

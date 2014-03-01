/*
 * Copyright 2014 Jonas Luethke
 */


package de.tuhh.luethke.okde.Exceptions;

public class EmptyDistributionException extends Exception {

	private static final long serialVersionUID = 2318130121620973423L;

	public EmptyDistributionException() {
	}

	public EmptyDistributionException(String msg) {
		super(msg);
	}
}

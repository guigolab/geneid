package gphase.tools;

import java.util.EmptyStackException;

public class Stack extends java.util.Stack {
	public boolean isEmpty() {
		return elementCount== 0;
	}
	
	public Object top() throws EmptyStackException {
		if (isEmpty())
			throw new EmptyStackException();
		return elementData[elementCount- 1];
	}
}

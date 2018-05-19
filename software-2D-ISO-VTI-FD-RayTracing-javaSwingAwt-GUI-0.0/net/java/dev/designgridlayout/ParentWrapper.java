//  Copyright 2005-2011 Jason Aaron Osgood, Jean-Francois Poilpret
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

package net.java.dev.designgridlayout;

import java.awt.Component;
import java.awt.Container;

import javax.swing.JComponent;

final class ParentWrapper<T extends Container>
{	
	ParentWrapper(T parent)
	{
		_parent = parent;
	}
	
	void checkParent(Container parent)
	{
		if (parent != _parent)
		{
			throw new IllegalArgumentException(
				"Must use layout instance with original parent container");
		}
	}
	
	void add(Component child)
	{
		try
		{
			_addChild = true;
			_parent.add(child);
		}
		finally
		{
			_addChild = false;
		}
	}
	
	void checkAddedComponent(JComponent component)
	{
		Container parent = component;
		while (parent != null)
		{
			if (parent == _parent)
			{
				throw new IllegalArgumentException("Do not add the same component twice");
			}
			parent = parent.getParent();
		}
	}

	void checkAdd()
	{
		if (!_addChild)
		{
			//TODO better message
			throw new IllegalArgumentException("Do not use this method");
		}
	}
	
	T parent()
	{
		return _parent;
	}
	
	final private T _parent;
	private boolean _addChild = false;
}

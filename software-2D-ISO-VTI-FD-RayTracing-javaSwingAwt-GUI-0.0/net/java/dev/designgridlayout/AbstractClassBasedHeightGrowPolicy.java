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

abstract class AbstractClassBasedHeightGrowPolicy<T extends Component>
	implements ClassBasedHeightGrowPolicy
{
	protected AbstractClassBasedHeightGrowPolicy(Class<T> componentClass)
	{
		_componentClass = componentClass;
	}

	@Override final public Class<? extends Component> getComponentClass()
	{
		return _componentClass;
	}

	@Override final public boolean canGrowHeight(Component component)
	{
		return componentCanGrowHeight(_componentClass.cast(component));
	}
	
	@Override final public int computeExtraHeight(Component component, int extraHeight)
	{
		return componentComputeExtraHeight(
			_componentClass.cast(component), extraHeight);
	}

	// Should be overridden if T components might not always be growable in height
	protected boolean componentCanGrowHeight(T component)
	{
		return true;
	}

	// Should be overridden if T components have special requirements in height
	protected int componentComputeExtraHeight(T component, int extraHeight)
	{
		return extraHeight;
	}

	private final Class<T> _componentClass;
}

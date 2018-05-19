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

class ContainerHeightGrowPolicy extends AbstractClassBasedHeightGrowPolicy<Container>
{
	public ContainerHeightGrowPolicy(HeightGrowPolicy defaultPolicy)
	{
		super(Container.class);
		_defaultPolicy = defaultPolicy;
	}

	@Override protected boolean componentCanGrowHeight(Container panel)
	{
		for (Component child: panel.getComponents())
		{
			if (_defaultPolicy.canGrowHeight(child))
			{
				return true;
			}
		}
		return false;
	}
	
	@Override protected int componentComputeExtraHeight(Container panel, int extraHeight)
	{
		int actualHeight = 0;
		for (Component child: panel.getComponents())
		{
			if (_defaultPolicy.canGrowHeight(child))
			{
				actualHeight = Math.max(actualHeight, 
					_defaultPolicy.computeExtraHeight(child, extraHeight));
			}
		}
		return actualHeight;
	}

	private final HeightGrowPolicy _defaultPolicy;
}

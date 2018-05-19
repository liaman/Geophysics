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

import javax.swing.JComponent;

import net.java.dev.designgridlayout.Componentizer.WidthPolicy;

// Used for all components added to a Componentizer
class ComponentizerItem extends BasicItem
{
	public ComponentizerItem(JComponent component, WidthPolicy widthPolicy)
	{
		super(component);
		_widthPolicy = widthPolicy;
	}

	@Override public int minimumWidth()
	{
		switch (_widthPolicy)
		{
			case PREF_FIXED:
			case PREF_AND_MORE:
			return preferredWidth();

			case MIN_TO_PREF:
			case MIN_AND_MORE:
			default:
			return super.minimumWidth();
		}
	}
	
	public WidthPolicy widthPolicy()
	{
		return _widthPolicy;
	}

	final private WidthPolicy _widthPolicy;
}

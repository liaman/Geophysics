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

abstract class AbstractRowItem extends BasicItem implements IRowItem
{
	protected AbstractRowItem(JComponent component)
	{
		super(component);
	}
	
	@Override public void hide()
	{
		if (isFirstSpanRow())
		{
			_isVisible =  component().isVisible();
			component().setVisible(false);
		}
	}

	@Override public void show()
	{
		if (isFirstSpanRow())
		{
			component().setVisible(_isVisible);
		}
	}

	private boolean _isVisible = true;
}

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

import javax.swing.JList;
import javax.swing.JScrollPane;

class JScrollPaneHeightGrowPolicy 
	extends AbstractClassBasedHeightGrowPolicy<JScrollPane>
{
	public JScrollPaneHeightGrowPolicy()
	{
		super(JScrollPane.class);
	}

	@Override protected int componentComputeExtraHeight(JScrollPane component, int extraHeight)
	{
		int unit = component.getVerticalScrollBar().getUnitIncrement(+1);
		// Fix for issue #28
		// TODO prepare a more extensible fix that can deal with any specific
		// component (only if needed: wait until other components require a fix)
		Component view = component.getViewport().getView();
		if (unit == 0 && view instanceof JList)
		{
			JList list = (JList) view;
			int visibleRows = list.getVisibleRowCount();
			if (visibleRows > 0)
			{
				unit = list.getPreferredScrollableViewportSize().height / visibleRows;
			}
		}
		// Make sure unit cannot be <= 0
		unit = Math.max(1, unit);
		// Return an integral number of units pixels
		return (extraHeight / unit) * unit;
	}
}

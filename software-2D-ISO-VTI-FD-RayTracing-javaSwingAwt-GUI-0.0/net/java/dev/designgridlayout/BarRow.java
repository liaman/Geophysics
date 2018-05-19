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

import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;

import javax.swing.JComponent;

final class BarRow extends AbstractRow implements IBarRow
{
	@Override public IBarRow add(JComponent child, Tag tag)
	{
		if (child != null)
		{
			checkAddedComponent(child);
			_items.add(new BarRowItem(child, (tag == null ? Tag.OTHER : tag)));
			parent().add(child);
		}
		return this;
	}

	@Override public IBarRow center(JComponent... children)
	{
		return add(Tag.OTHER, children);
	}

	@Override public IBarRow left(JComponent... children)
	{
		return add(Tag.LEFT, children);
	}

	@Override public IBarRow right(JComponent... children)
	{
		return add(Tag.RIGHT, children);
	}
	
	private IBarRow add(Tag tag, JComponent... children)
	{
		for (JComponent child: children)
		{
			add(child, tag);
		}
		return this;
	}
	
	@Override public IBarRow gap()
	{
		_items.add(null);
		return this;
	}

	@Override public IBarRow withOwnRowWidth()
	{
		_ownRowWidth = true;
		return this;
	}

	@Override List<BarRowItem> items()
	{
		if (_leftItems == null)
		{
			// First of all split into 3 lists: left, center and right
			_leftItems = PlatformHelper.extractLeftItems(_items);
			_centerItems = PlatformHelper.extractCenterItems(_items);
			_rightItems = PlatformHelper.extractRightItems(_items);

			// Then rebuild _items to only include all but non-null items
			_items.clear();
			_items.addAll(_leftItems);
			_items.addAll(_centerItems);
			_items.addAll(_rightItems);
			// Calculate total number of extra gaps
			_numUnrelatedGaps = numGapsBetweenParts();
			Iterator<BarRowItem> i = _items.iterator();
			while (i.hasNext())
			{
				if (i.next() == null)
				{
					i.remove();
					_numUnrelatedGaps++;
				}
			}
		}
		return _items;
	}

	@Override int totalNonGridWidth(int hgap, int unrelhgap)
	{
		int leftWidth = computePartWidth(_leftItems, hgap, unrelhgap);
		int rightWidth = computePartWidth(_rightItems, hgap, unrelhgap);
		int centerWidth = computePartWidth(_centerItems, hgap, unrelhgap);
		int sidesWidth = Math.max(leftWidth, rightWidth);
		if (centerWidth != 0 && (leftWidth != 0 || rightWidth != 0))
		{
			sidesWidth *= 2;
		}
		int totalWidth = centerWidth + sidesWidth + unrelhgap * numGapsBetweenParts();
		return totalWidth;
	}

	@Override int componentNonGridWidth()
	{
		return (_ownRowWidth ? 0 : actualComponentNonGridWidth());
	}
	
	private int actualComponentNonGridWidth()
	{
		return ComponentHelper.maxValues(_items, PrefWidthExtractor.INSTANCE);
	}

	@Override void forceComponentNonGridWidth(int width)
	{
		_compWidth = ((width > 0 && !_ownRowWidth) ? width : actualComponentNonGridWidth());
	}

	//CSOFF: ParameterNumber
	@Override int layoutRow(LayoutHelper helper, int left, int hgap, int gridgap, 
		int unrelhgap, int rowWidth, int gridsWidth, List<Integer> labelsWidth)
	{
		// Layout each part of the row individually
		int x = left;
		int actualHeight = layoutOnePart(helper, x, hgap, unrelhgap, _leftItems);

		x = left + (rowWidth - computePartWidth(_centerItems, hgap, unrelhgap)) / 2;
		actualHeight = Math.max(
			actualHeight, layoutOnePart(helper, x, hgap, unrelhgap, _centerItems));

		x = left + rowWidth - computePartWidth(_rightItems, hgap, unrelhgap);
		actualHeight = Math.max(
			actualHeight, layoutOnePart(helper, x, hgap, unrelhgap, _rightItems));

		return actualHeight;
	}
	//CSON: ParameterNumber
	
	private int computePartWidth(List<BarRowItem> items, int hgap, int unrelhgap)
	{
		if (!items.isEmpty())
		{
			int numUnrelGaps = Collections.frequency(items, null);
			int numComponents = items.size() - numUnrelGaps;
			int numGaps = numComponents - numUnrelGaps - 1;
			int width = numComponents * _compWidth + numGaps * hgap + numUnrelGaps * unrelhgap;
			return width;
		}
		else
		{
			return 0;
		}
	}
	
	private int numGapsBetweenParts()
	{
		// Calculate total number of extra gaps
		int numParts =	(_leftItems.isEmpty() ? 0 : 1) +
						(_centerItems.isEmpty() ? 0 : 1) +
						(_rightItems.isEmpty() ? 0 : 1);
		return (numParts > 0 ? numParts - 1 : 0);
	}

	private int layoutOnePart(
		LayoutHelper helper, int xOrigin, int hgap, int unrelhgap, List<BarRowItem> items)
	{
		int x = xOrigin;
		int actualHeight = 0;
		for (BarRowItem item: items)
		{
			if (item != null)
			{
				actualHeight = Math.max(actualHeight, helper.setSizeLocation(
					item.component(), x, _compWidth, height(), baseline()));
				x += _compWidth + hgap;
			}
			else
			{
				x += unrelhgap - hgap;
			}
		}
		return actualHeight;
	}
	
	private final List<BarRowItem> _items = new ArrayList<BarRowItem>();
	private List<BarRowItem> _leftItems = null;
	private List<BarRowItem> _centerItems = null;
	private List<BarRowItem> _rightItems = null;
	private int _numUnrelatedGaps = 0;
	private boolean _ownRowWidth = false;
	private int _compWidth = 0;
}

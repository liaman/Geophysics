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

import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

final class RowIterator implements Iterator<AbstractRow>
{
	static public Iterable<AbstractRow> each(final List<AbstractRow> rows)
	{
		return new Iterable<AbstractRow>()
		{
			@Override public Iterator<AbstractRow> iterator()
			{
				return new RowIterator(rows);
			}
		};
	}

	private RowIterator(List<AbstractRow> rows)
	{
		_rows = rows;
	}
	
	@Override public AbstractRow next()
	{
		int index = findNext();
		if (index != -1)
		{
			_index = index + 1;
			return _rows.get(index);
		}
		else
		{
			throw new NoSuchElementException();
		}
	}
	
	@Override public boolean hasNext()
	{
		return findNext() != -1;
	}
	
	@Override public void remove()
	{
		throw new UnsupportedOperationException();
	}
	
	private int findNext()
	{
		for (int i = _index; i < _rows.size(); i++)
		{
			AbstractRow row = _rows.get(i);
			if (!row.isEmpty())
			{
				return i;
			}
		}
		return -1;
	}
	
	private final List<AbstractRow> _rows;
	private int _index = 0;
}
